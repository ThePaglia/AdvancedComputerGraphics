#version 420

layout(location = 0) out vec4 fragmentColor;

uniform float uTime;
uniform vec2 uResolution;
uniform sampler2D uNoiseTexture;

// Camera
uniform vec3 uCameraPos;
uniform vec3 uCameraRight;
uniform vec3 uCameraUp;
uniform vec3 uCameraDir;

// Light
uniform vec3 pointLightColor;
uniform vec3 lightPosition;
uniform float pointLightIntensityMultiplier;

// Cloud
const vec3 ambientColor = vec3(0.60, 0.60, 0.75);
const float ambientIntensity = 0.9f;
uniform float cloudTime = 1.0f;
uniform float cloudShadowIntensity = 4f;
uniform float cloudShadowCutoff = 0.5f;


// Raymarching
#define MAX_STEPS 500
#define MAX_MARCH_DISTANCE 1000
const float MARCH_SIZE = 0.08;
const float ATMOSPHERE_MARCH_SIZE = 100;
uniform sampler2D uBlueNoise;

// Sampling
uniform float samplingIncreaseFactor;
uniform float samplingIncreaseDepth;
uniform float samplingFalloffDistance;

// Scene parameters
const float cloudHeight = 10f;
const float cloudBoxWidth = 200;
const float atmosphereFalloffDepth = 100;

// Planet parameters
uniform vec3 planetOrigin = vec3(0, 0, 0);
uniform float planetRadius = 10f;
uniform float cloudlessDepth = 0.5f;
uniform float cloudDepth = 1f;
uniform float cloudScale = 1.0f;
uniform float atmosphereDepth = 3f;
uniform float cloudStepMin = 0.1f;
uniform float cloudStepMax = 0.8f;
uniform float atmosphereDensityFalloff = 2f;
uniform vec3 atmosphereScatteringCoefficients = vec3(0, 0, 0);

// Precalculated constants
const float atmosphereRadius = planetRadius + atmosphereDepth;
vec3 sunDirection = normalize(lightPosition);

float sdSphere(vec3 p, float radius) {
	return length(p) - radius;
}

float noise(vec3 x) {
	vec3 p = floor(x);
	vec3 f = fract(x);
	f = f * f * (3.0 - 2.0 * f);

	vec2 uv = (p.xy + vec2(37.0, 239.0) * p.z) + f.xy;
	vec2 tex = textureLod(uNoiseTexture, (uv + 0.5) / 256.0, 0.0).yx;

	return mix(tex.x, tex.y, f.z) * 2.0 - 1.0;
}

// Fractal Brownian Motion
float fbm(vec3 p) {
	vec3 q = p + cloudTime * vec3(1.0, -0.2, -1.0);
	float g = noise(q);

	float f = 0.0;
	float scale = 0.5;
	float factor = 2.02;

	for(int i = 0; i < 6; i++) {
		f += scale * noise(q);
		q *= factor;
		factor += 0.21;
		scale *= 0.5;
	}

	return f;
}

float evaluateDensityAt(vec3 p) {
    float f = fbm(p * 3);

    //float cloudBelt = min(-p.y - (-cloudHeight) + cloudDepth / 2, p.y + (-cloudHeight) + cloudDepth / 2) + f;
    
    float outerCloudShell = sdSphere(p - planetOrigin, planetRadius + cloudlessDepth + cloudDepth);
    float innerCloudShell = -sdSphere(p - planetOrigin, planetRadius + cloudlessDepth); // cloudDepth = shell thickness

    float shell = max(outerCloudShell, innerCloudShell); // Hollow region
    float shellDensity = -shell;

    // Large-scale holes
    // NOTE! Something fishy happens when cloudScale goes beyond 1, you get these black clumps of darkness which I believe come from either how we scale the clouds, or how the smoothstep works
    float largeHoleNoise = fbm(p * cloudScale); // lower frequency = larger features
    float holeMask = smoothstep(cloudStepMin, cloudStepMax, largeHoleNoise); // controls size/sharpness of holes

    return holeMask * shellDensity;
}

float sdf(vec3 p) {
    float sphere = sdSphere(p - planetOrigin, planetRadius);
    return sphere;
}

vec3 calculateNormal(in vec3 p)
{
    const vec3 small_step = vec3(0.001, 0.0, 0.0);

    float gradient_x = sdf(p + small_step.xyy) - sdf(p - small_step.xyy);
    float gradient_y = sdf(p + small_step.yxy) - sdf(p - small_step.yxy);
    float gradient_z = sdf(p + small_step.yyx) - sdf(p - small_step.yyx);

    vec3 normal = vec3(gradient_x, gradient_y, gradient_z);

    return normalize(normal);
}


// Stolen from https://www.shadertoy.com/view/Ml3Gz8
float smoothmin(float a, float b, float k) {
    
    // Compute the difference between the two values.
    // This is used to interpolate both values inside the range (-k, k).
    // Smaller ranges give a better approximation of the min function.
    float h = a - b;
    
    // The interval [-k, k] is mapped to [0, 1],
    // and clamping takes place only after this transformation.
    
    // Map [-k, k] to [0, 1] and clamp if outside the latter.
    h = clamp(0.5 + 0.5*h/k, 0.0, 1.0);    
    
    // Linearly interpolate the input values using h inside (0, 1).
    // The second term ensures continuous derivatives at the boundaries of [0,1],
    // but this is not completely obvious! See my blog post for details.
    return mix(a, b, h) - k*h*(1.0-h);    
}

bool intersectBox(vec3 ro, vec3 rd, vec3 boxMin, vec3 boxMax, out float tEnter, out float tExit) {
    vec3 invDir = 1.0 / rd;
    vec3 t0 = (boxMin - ro) * invDir;
    vec3 t1 = (boxMax - ro) * invDir;

    vec3 tMin = min(t0, t1);
    vec3 tMax = max(t0, t1);

    tEnter = max(max(tMin.x, tMin.y), tMin.z);
    tExit  = min(min(tMax.x, tMax.y), tMax.z);

    return tEnter <= tExit;
}

bool intersectSphere(vec3 ro, vec3 rd, vec3 sc, float sr, out float tEnter, out float tExit) {
    vec3 oc = ro + sc;
    float b = dot(oc, rd);
    float c = dot(oc, oc) - sr * sr;
    float h = b * b - c;

    if (h < 0.0) return false;  // No intersection

    h = sqrt(h);
    tEnter = -b - h;
    tExit  = -b + h;
    return true;
}

bool raymarchSDF(vec3 rayOrigin, vec3 rayDirection, int maxSteps, float maxDepth, out float depth, out vec3 point) {
    for(int i = 0; i < maxSteps; i++) {
        point = rayOrigin + rayDirection * depth;
        float distanceToWorld = sdf(point);

        if(distanceToWorld <= 0.01f) {
            return true;    
        }

        depth += distanceToWorld;

        if(depth > maxDepth) {
            return false;
        }
    }

    return false;
}

// Atmospheric density has an exponential falloff the higher we get from the surface
float atmosphereDensityAtPoint(vec3 p) {
    float heightAboveSurface = sdSphere(p - planetOrigin, planetRadius);
    // normalized height from 1 at the edge of the atmosphere shell to 0 at the surface
    float normalizedHeight = heightAboveSurface / atmosphereDepth;
    // Multiplying by (1 - normalizedHeight) is a bit dirty but ensures that the density is 0 at the edge of the atmosphere shell
    float density = exp(-normalizedHeight * atmosphereDensityFalloff) * (1 - normalizedHeight);
    return density;
}

float numOpticalDepthPoints = 10;
float atmosphereOpticalDepth(vec3 rayOrigin, vec3 rayDirection, float rayLength) {
    float stepSize = rayLength / (numOpticalDepthPoints - 1);
    float opticalDepth = 0.0f;
    
    float marchDepth = 0.0;

    // March through the atmosphere along the ray length and accumulate density
    for(int i = 0; i < numOpticalDepthPoints; i++) {
        vec3 p = rayOrigin + rayDirection * marchDepth;

        float density = atmosphereDensityAtPoint(p);
        opticalDepth += density * stepSize;
        marchDepth += stepSize;
    }

    return opticalDepth;
}

float numInscatteringPoints = 10;
vec4 calculateAtmosphereLight(vec3 rayOrigin, vec3 rayDirection, float rayLength, vec4 originalColor) {
    float stepSize = rayLength / (numInscatteringPoints - 1);
    vec3 inScatteredLight = vec3(0.0);
    float viewRayOpticalDepth = 0.0;

    float marchDepth = 0.0;

    for(int i = 0; i < numInscatteringPoints; i++) {
        vec3 p = rayOrigin + rayDirection * marchDepth;
        // The distance of the sun ray, i.e. the ray that travels from p through the atmosphere towards the sun, is the same as tExit
        float tEnter, sunRayLength;
        // We are guaranteed to intersect the atmosphere sphere as we start the ray inside of it, we can therefore use this function to get sunRayLength
        intersectSphere(p, sunDirection, planetOrigin, atmosphereRadius, tEnter, sunRayLength);
        
        // The accumulated amount of light (density) from the point towards the sun
        float sunRayOpticalDepth = atmosphereOpticalDepth(p, sunDirection, sunRayLength);
        
        // The accumulated amount of light (density) from the point back towards the camera
        // NOTE: this variable creates a somewhat noticeable ring of darkness right around the planet as this is where a ray travels the furthest through the atmosphere, 
        // this behaviour seems to me to be physically correct but the result is somewhat strange
        // It also seems to lead to a pronounciation of the atmosphere's color bands
        viewRayOpticalDepth = atmosphereOpticalDepth(p, -rayDirection, marchDepth);

        // The light that reaches the camera has an exponential falloff
        vec3 transmittance = exp(-(sunRayOpticalDepth + viewRayOpticalDepth) * atmosphereScatteringCoefficients);
        
        // Sample the density at this inscatter point
        float density = atmosphereDensityAtPoint(p);

        // Accumulate the density multiplied by the transmittance at this point, i.e. the amount of light that reaches this point and is transmitted towards the camera 
        inScatteredLight += density * transmittance * atmosphereScatteringCoefficients * stepSize;

        marchDepth += stepSize;
    }

    // TODO: Figure out how to properly simulate how much of the light from the planet's surface that is scattered away from the camera
    float originalColorTransmittance = 1;

    return originalColor * originalColorTransmittance + vec4(inScatteredLight, 0);
}

vec4 raymarch(vec3 rayOrigin, vec3 rayDirection, vec3 cameraForward, float offset) {
    float rayDotCam = dot(rayDirection, cameraForward);
    float rayDotCloudPlane = dot(rayDirection, normalize(rayDirection * vec3(1, 0, 1)));
    float distCamCloudPlane = abs(rayOrigin.y - cloudHeight);

    vec4 opaqueRes = vec4(0.0);

    float opaqueDepth = 0;
    vec3 opaquePoint = vec3(0);
    if(raymarchSDF(rayOrigin, rayDirection, MAX_STEPS, MAX_MARCH_DISTANCE, opaqueDepth, opaquePoint)) {
        vec3 normal = calculateNormal(opaquePoint);
        float diffuseIntensity = clamp(dot(sunDirection, normal), 0, 1);

        float cloudDensityAbove = 0;

        float tEnterInner, tExitInner;
        int numShadowSteps = 8;
        if(diffuseIntensity > 0 && intersectSphere(opaquePoint, sunDirection, planetOrigin, planetRadius + cloudlessDepth, tEnterInner, tExitInner)) {
            float tEnterOuter, tExitOuter;
            if(intersectSphere(opaquePoint, sunDirection, planetOrigin, planetRadius + cloudlessDepth + cloudDepth, tEnterOuter, tExitOuter)) {
                float rayLength = tExitOuter - tExitInner;
                float stepSize = rayLength / (numShadowSteps - 1);
                for(int i = 0; i < numShadowSteps; i++) {
                    // You could add the offset here to make shadow artifacting somewhat less noticeable
                    vec3 p = opaquePoint + sunDirection * (tExitInner + i * stepSize);

                    float density = evaluateDensityAt(p);
                    cloudDensityAbove += max(density, 0) * stepSize;
                }
            }
        }

        cloudDensityAbove *= cloudShadowIntensity;

        opaqueRes = mix(vec4(0, 0, 0, 1), vec4(0.01f, 0.2f, 1.0f, 1.0f), diffuseIntensity);
        opaqueRes.rgb *= pointLightColor * pointLightIntensityMultiplier;
        opaqueRes.rgb *= clamp(1 - cloudDensityAbove, cloudShadowCutoff, 1.0);
        opaqueRes = pow(opaqueRes, vec4(1 / 2.2f)); // Gamma correction
    }

    // Right now the atmosphere is rendered on top of the clouds as their depth cannot be used to cut off the distanceThroughAtmosphere
    // Perhaps some clever compositing can save us here?
    vec4 atmosphereLight = vec4(0.0f);
    float atmosphereRadius = planetRadius + atmosphereDepth;
    float tEnterAtmosphere, tExitAtmosphere;
    if(intersectSphere(rayOrigin, rayDirection, planetOrigin, atmosphereRadius, tEnterAtmosphere, tExitAtmosphere) && tExitAtmosphere > 0) {
        vec3 p = rayOrigin + rayDirection * max(tEnterAtmosphere, 0);
        float distanceThroughAtmosphere = min(tExitAtmosphere, opaqueDepth) - max(tEnterAtmosphere, 0);
        atmosphereLight = calculateAtmosphereLight(p, rayDirection, distanceThroughAtmosphere, opaqueRes);
    }

    vec3 closestPointToSunOnAtmosphereShell = planetOrigin + sunDirection * atmosphereRadius;

    float volumetricDepth = 0;
    vec4 volumetricRes = vec4(0.0);
    float tEnterClouds, tExitClouds;
    if(intersectSphere(rayOrigin, rayDirection, planetOrigin, planetRadius + cloudlessDepth + cloudDepth, tEnterClouds, tExitClouds) && tExitClouds > 0) {
        float tEnterInnerSphere, tExitInnerSphere;
        if(!intersectSphere(rayOrigin, rayDirection, planetOrigin, planetRadius + cloudlessDepth, tEnterInnerSphere, tExitInnerSphere))
            tEnterInnerSphere = MAX_MARCH_DISTANCE;
        float startDepth = max(tEnterClouds, 0.0);
        float camAlignedPlane = ceil((startDepth * rayDotCam) / MARCH_SIZE) * MARCH_SIZE;
        volumetricDepth = camAlignedPlane / rayDotCam;  // World-space depth, but camera-aligned
        startDepth = volumetricDepth;
    
        float depthTraveledThroughMedium = 0;
        bool hasJumpedForward = false;
        for(int i = 0; i < MAX_STEPS; i++) {

            vec3 p = rayOrigin + rayDirection * volumetricDepth;

            float density = evaluateDensityAt(p);

            if (density > 0.0) {
                float shadowMultiplier = 1;
                float ent, ext;

                // Hacky way of determining if a cloud is in shadow, should be updated to be physically correct in the future
                // Uses a similar technique to lambertian diffuse lighting, but after the top half of the planet (top half is fully lit)
                // The bottom half has a gradient from 1 at the middle to 0 at the bottom, but the gradient can be adjusted:
                // lightingFalloff adjusts the cutoff point of the gradient of the bottom half of the planet, for example:
                // 1.0 -> gradient extends the whole bottom half
                // 0.5 -> gradient extends half of the bottom half
                // 0.1 -> gradient extends one 10th of the bottom half
                float lightingFalloff = 0.5;
                // TODO: calculate lightingFalloff from the difference of the cloudRadius and the atmosphereRadius
                float diffuseIntensity = clamp(dot(sunDirection, normalize(p - planetOrigin)) * (1 / lightingFalloff) + 1, 0, 1);

                // Scale cloud lighting by how far away it is from the closest point on the atmosphere shell, 
                // hacky and not (even close to) a physically correct way of achieving darker clouds at the far end of the planet
                shadowMultiplier *= diffuseIntensity;

                float diffuse = clamp((density - evaluateDensityAt(p + 0.3 * sunDirection)) / 0.3, 0.0, 1.0);
                // TODO: Make cloud color reflect the incoming sunlight's color, i.e. a nice orange/red at glancing angles
                vec3 lin = ambientColor * ambientIntensity + pointLightIntensityMultiplier * pointLightColor * diffuse;
                vec4 color = vec4(mix(vec3(1.0), vec3(0.0), density), density);
                color.rgb *= lin;
                color.rgb *= color.a;
                color.rgb *= shadowMultiplier;
                // TODO: calculate the scattering of the cloud's light as it travels towards the camera
                //color *= exp(-atmosphereOpticalDepth(p, -rayDirection, volumetricDepth - tEnterAtmosphere));
                volumetricRes += color * (1.0 - volumetricRes.a);
                // This doesn't seem to be enough to completely get rid of the outline of the planet showing through clouds, as the atmosphere stepping size is shorter when stepping towards the edge of the planet compared to when stepping just beyond it
                // Perhaps it would be possible to also march the atmosphere by a constant amount? Perhaps this would mitigate the issue?
                if (volumetricRes.a >= 0.999) {
                    break;
                }
            }

            depthTraveledThroughMedium += MARCH_SIZE;

            // Samples should be taken at higher frequencies when closer to the camera, i.e. when the depth is low
            float depthFactor = clamp(max(volumetricDepth - samplingIncreaseDepth, 0) / MAX_MARCH_DISTANCE, 0, 1) * clamp(samplingFalloffDistance / max(distCamCloudPlane, 1), 0, 1);
            volumetricDepth = startDepth + depthTraveledThroughMedium * (samplingIncreaseFactor * depthFactor  + 1);

            // TODO: Figure out why we don't get a performance increase when jumping forward through the empty space between the planet and the start of the cloud layer (which is what this part of the code is meant to do)
            // NOTE: I realized it's because we don't update startDepth as well! -> PLEASE FIX
            if(!hasJumpedForward && volumetricDepth > tEnterInnerSphere) {
                volumetricDepth = tExitInnerSphere;
                hasJumpedForward = true;
            }

            if(volumetricDepth >= opaqueDepth || volumetricDepth >= tExitClouds || depthTraveledThroughMedium > MAX_MARCH_DISTANCE)
                break;
            
        }
    }

    vec4 res = volumetricRes + atmosphereLight * (1 - volumetricRes.a);

    return res;
}

void main() {
    sunDirection = normalize(lightPosition);
	mat3 uCameraMatrix = transpose(mat3(uCameraRight, uCameraUp, uCameraDir));
	vec2 uv = gl_FragCoord.xy / uResolution.xy;
	uv -= 0.5;
	uv.x *= uResolution.x / uResolution.y;
	
	// Ray Origin - camera
	vec3 rayOrigin = uCameraPos;
	// Ray Direction
	vec3 rayDirection = normalize(vec3(uv, -1.0) * uCameraMatrix);

	vec3 color = vec3(0.0);
	float blueNoise = texture2D(uBlueNoise, gl_FragCoord.xy / 100).r;
	float offset = fract(blueNoise);
	vec4 res = raymarch(rayOrigin, rayDirection, uCameraDir, offset * 0.1f);
	color = res.rgb;

	fragmentColor = vec4(color, 1.0);
}
