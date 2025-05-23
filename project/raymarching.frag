#version 420

in vec2 vUV;
layout(location = 0) out vec4 fragmentColor;

uniform float uTime;
uniform vec2 uResolution;
uniform sampler2D uSceneColor; // Result from rasterization step
uniform sampler2D uSceneDepth; // Result from rasterization step
uniform sampler2D uNoiseTexture;

// These are values retrieved from the previous rasterized step
vec3 sceneColor;
float sceneDepth;
vec3 scenePoint;

// Camera
uniform vec3 uCameraPos;
uniform vec3 uCameraRight;
uniform vec3 uCameraUp;
uniform vec3 uCameraDir;
uniform mat4 uViewProjectionMatrix;

// Light
uniform vec3 directionalLightColor;
uniform vec3 lightPosition;
uniform float directionalLightIntensityMultiplier;

// Cloud
const vec3 ambientColor = vec3(0.60, 0.60, 0.75);
const float ambientIntensity = 0.9f;
uniform float cloudTime = 1.0f;
uniform float cloudShadowIntensity = 4f;
uniform float cloudShadowCutoff = 0.5f;
uniform float cloudNoiseUVScale = 128.0f;
uniform float cloudNoiseAmount = 0.1f;

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
uniform float cloudLightingFalloff = 0.5f;
uniform float atmosphereDensityFalloff = 2f;
uniform vec3 atmosphereScatteringCoefficients = vec3(0, 0, 0);
uniform float atmosphereDensityAtSeaLevel = 0.5f;

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
    vec3 oc = ro - sc;
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
    float density = exp(-normalizedHeight * atmosphereDensityFalloff) * (1 - normalizedHeight) * atmosphereDensityAtSeaLevel;
    return density;
}

float numOpticalDepthPoints = 8;
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

float numInscatteringPoints = 8;
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

    vec4 opaqueRes = vec4(sceneColor, 1);

    // Calculate the pixel contribution of opaque geometry
    vec3 opaquePoint = vec3(scenePoint);
    float opaqueDepth = length(rayOrigin - opaquePoint);
    if(sceneColor != vec3(0)) {
        float cloudDensityAbove = 0;

        // Calculate cloud density above (i.e. in the sun's direction) this point
        float tEnterInner, tExitInner;
        int numShadowSteps = 8;
        if(intersectSphere(opaquePoint, sunDirection, planetOrigin, planetRadius + cloudlessDepth, tEnterInner, tExitInner)) {
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

        opaqueRes.rgb *= clamp(1 - cloudDensityAbove, cloudShadowCutoff, 1.0);
    }

    // Calculate the atmosphere lighting contribution
    vec4 atmosphereLight = vec4(0.0f);
    float atmosphereRadius = planetRadius + atmosphereDepth;
    float tEnterAtmosphere, tExitAtmosphere;
    if(intersectSphere(rayOrigin, rayDirection, planetOrigin, atmosphereRadius, tEnterAtmosphere, tExitAtmosphere) && tExitAtmosphere > 0) {
        vec3 p = rayOrigin + rayDirection * max(tEnterAtmosphere, 0);
        float distanceThroughAtmosphere = min(tExitAtmosphere, opaqueDepth) - max(tEnterAtmosphere, 0);
        atmosphereLight = calculateAtmosphereLight(p, rayDirection, distanceThroughAtmosphere, opaqueRes);
    }
    
    // The depth of the current ray marched step through the clouds
    float volumetricDepth = 0;
    vec4 volumetricRes = vec4(0.0);

    // Accumulate cloud color (volumetricRes) if we hit the cloud shell
    float tEnterClouds, tExitClouds;
    if(intersectSphere(rayOrigin, rayDirection, planetOrigin, planetRadius + cloudlessDepth + cloudDepth, tEnterClouds, tExitClouds) && tExitClouds > 0) {
        float startDepth = max(tEnterClouds, 0.0);
        
        // We need to keep track of whether or not we have skipped the inner cloud shell, i.e. entered a cloud region and then exited it and jumped forward to the next cloud region on the other side of the cloud shell
        bool hasSkippedInnerCloudShell = false;

        // Check if we hit the space in between the clouds and the planet (the "inner sphere"). If so, we can skip traversing it as there are no clouds there
        float tEnterInnerSphere, tExitInnerSphere;
        if(intersectSphere(rayOrigin, rayDirection, planetOrigin, planetRadius + cloudlessDepth, tEnterInnerSphere, tExitInnerSphere) && tExitInnerSphere > 0) {
            // If we are inside the space between the clouds and the planet (i.e. if tEnterInnerSphere < 0, i.e. we enter the sphere behind us), jump forward to just before we enter the clouds
            startDepth = tEnterInnerSphere <= 0 ? max(tExitInnerSphere, 0.0) - MARCH_SIZE : startDepth;
            // If we are inside the space between the clouds and the planet, set hasSkippedInnerCloudShell to true, this ensures that we won't jump forward later again (as we have already done it on the last row)
            hasSkippedInnerCloudShell = tEnterInnerSphere <= 0;
        } else {
            // If we did not hit the space between the clouds and the planet, set hasSkippedInnerCloudShell to true, this ensures that we won't jump forward later
            hasSkippedInnerCloudShell = true;
        }

        // Set the start depth to a camera aligned plane. This gets rid of "hitbox" artifacts arising from hitting the cloud sphere at non-discretisized distances from the camera
        float camAlignedPlane = floor((startDepth * rayDotCam) / MARCH_SIZE) * MARCH_SIZE;
        volumetricDepth = camAlignedPlane / rayDotCam - offset * cloudNoiseAmount * (tEnterClouds < 0 ? -1 : 1);  // World-space depth, but camera-aligned
        startDepth = volumetricDepth;

        // This is used to keep track of the distance traveled through the clouds
        float depthTraveledThroughMedium = 0;

        // Raymarch through the cloud shell and accumulate the result
        for(int i = 0; i < MAX_STEPS; i++) {
            // Get the point at the current depth
            vec3 p = rayOrigin + rayDirection * volumetricDepth;

            // Get the density at the current point
            float density = evaluateDensityAt(p);

            // Accumulate the result if the density is above 0.0
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
                float lightingFalloff = cloudLightingFalloff;
                // TODO: calculate lightingFalloff from the difference of the cloudRadius and the atmosphereRadius
                float diffuseIntensity = clamp(dot(sunDirection, normalize(p - planetOrigin)) * (1 / lightingFalloff) + 1, 0, 1);

                // Scale cloud lighting by how far away it is from the closest point on the atmosphere shell, 
                // hacky and not (even close to) a physically correct way of achieving darker clouds at the far end of the planet
                shadowMultiplier *= diffuseIntensity;

                // TODO: perhaps it is possible to increment the opticalDepth using atmosphereDensityAtPoint as we step through the atmosphere? Instead of recalculating it every frame (which the following line does)
                float viewRayOpticalDepth = atmosphereOpticalDepth(p, -rayDirection, volumetricDepth - max(tEnterAtmosphere, 0));

                float diffuse = clamp((density - evaluateDensityAt(p + 0.3 * sunDirection)) / 0.3, 0.0, 1.0);
                // TODO: Make cloud color reflect the incoming sunlight's color, i.e. a nice orange/red at glancing angles
                vec3 lin = ambientColor * ambientIntensity + directionalLightIntensityMultiplier * directionalLightColor * diffuse;
                vec4 color = vec4(mix(vec3(1.0), vec3(0.0), density), density);
                color.rgb *= lin;
                color.rgb *= color.a;
                color.rgb *= shadowMultiplier;
                color *= exp(-viewRayOpticalDepth); // Not sure if this is the best way of multiplying the contribution of the viewRayOpticalDepth as it affects the transparency
                volumetricRes += color * (1.0 - volumetricRes.a);

                // We can immediately break out of the loop if the transparency is greater than this treshold, the reasoning is that any further steps would contribute an insignificant amount to the pixel color
                if (volumetricRes.a >= 0.999) {
                    break;
                }
            }

            // Move forward one step through the medium
            depthTraveledThroughMedium += MARCH_SIZE;
            volumetricDepth = startDepth + depthTraveledThroughMedium;

            // Jump forward through the empty space between the planet and the start of the cloud layer if we enter the "inner shell" (as described before)
            // NOTE: that jumping forward still leads to some artifacting around the edges of the inner sphere, at glancing angles you can kind of see where the shell starts. Probably has something to do with how we jump forward
            if(!hasSkippedInnerCloudShell && volumetricDepth >= tEnterInnerSphere) {
                startDepth += tExitInnerSphere - max(tEnterInnerSphere, 0.0) - MARCH_SIZE;

                // We must align the new starting point to a discretisized value from the camera plane
                float camAlignedPlane = floor((startDepth * rayDotCam) / MARCH_SIZE) * MARCH_SIZE;
                volumetricDepth = camAlignedPlane / rayDotCam + offset * cloudNoiseAmount * (tEnterClouds < 0 ? -1 : 1);  // World-space depth, but camera-aligned
                startDepth = volumetricDepth;
                hasSkippedInnerCloudShell = true;
            }

            // We can break out of the loop if we hit the opaqueDepth or exit the outer edge of the cloud shell
            if(volumetricDepth >= opaqueDepth || volumetricDepth >= tExitClouds)
                break;
        }
    }

    // Blend together the atmospheric result with the volumetric one, depending on volumetric alpha
    vec4 res = volumetricRes + atmosphereLight * (1 - volumetricRes.a);

    return res;
}

void main() {
    vec2 uv = gl_FragCoord.xy / uResolution.xy;
    sceneColor = texture(uSceneColor, uv).rgb;
    float depth = texture(uSceneDepth, uv).r;

    sunDirection = normalize(lightPosition);

    // Next, we must construct the correct ray direction from the same view projection matrix that the rasterized geometry uses

    // 1. Normalized Device Coordinates (NDC)
    vec2 normalizedDeviceCoordinates = uv * 2.0 - 1.0;

    // 2. Create clip-space position
    float z = depth * 2.0 - 1.0;
    vec4 clipSpacePos = vec4(normalizedDeviceCoordinates, z, 1.0);

    // 3. Transform by inverse of viewProjectionMatrix to get world-space position
    vec4 worldPoint = inverse(uViewProjectionMatrix) * clipSpacePos;
    worldPoint /= worldPoint.w;
    // Store the scene point as a vec3, as we no longer need information about w
    scenePoint = vec3(worldPoint);
    sceneDepth = length(scenePoint - uCameraPos);

    // 4. Ray origin is camera position
    vec3 rayOrigin = uCameraPos;

    // 5. Ray direction from camera to worldPos
    vec3 rayDirection = normalize(scenePoint.xyz - rayOrigin);

    // 6. Raymarching
    vec3 color = vec3(0.0);
    float blueNoise = texture2D(uNoiseTexture, gl_FragCoord.xy / cloudNoiseUVScale).r;
    float offset = fract(blueNoise);
    vec4 res = raymarch(rayOrigin, rayDirection, uCameraDir, offset);
    color = res.rgb;

    fragmentColor = vec4(color, 1.0);
}
