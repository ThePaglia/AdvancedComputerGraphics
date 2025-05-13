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
uniform vec3 planetOrigin = vec3(0, 0, -20);
uniform float planetRadius = 10f;
uniform float cloudlessDepth = 0.5f;
uniform float cloudDepth = 1f;
uniform float cloudScale = 1.0f;
uniform float cloudStepMin = 0.1f;
uniform float cloudStepMax = 0.8f;

struct cloud {
    vec3 position;
    float radius;
};

cloud[5] clouds = { { vec3(0, 0, 0), 2 }, { vec3(10, 0, 0), 5 }, { vec3(60, 0, 0), 2 }, { vec3(100, 0, 0), 8 }, { vec3(30, 0, 0), 1 } };

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

// Constructs a unified min from the list of clouds
float evaluateDensityAt(vec3 p) {
    float f = fbm(p * 3);

    float cloudBelt = min(-p.y - (-cloudHeight) + cloudDepth / 2, p.y + (-cloudHeight) + cloudDepth / 2) + f;
    
    float outer = sdSphere(p + planetOrigin, planetRadius + cloudlessDepth + cloudDepth);
    float inner = -sdSphere(p + planetOrigin, planetRadius + cloudlessDepth); // cloudDepth = shell thickness

    float shell = max(outer, inner); // Hollow region
    float shellDensity = -shell;

    // Large-scale holes
    float largeHoleNoise = fbm(p * cloudScale); // lower frequency = larger features
    float holeMask = smoothstep(cloudStepMin, cloudStepMax, largeHoleNoise); // controls size/sharpness of holes

    return shellDensity * holeMask;
}

float sdf(vec3 p) {
    float sphere = sdSphere(p + planetOrigin, planetRadius);
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

float remap(float value, float inMin, float inMax, float outMin, float outMax) {
    return outMin + (value - inMin) * (outMax - outMin) / (inMax - inMin);
}

float sampleAtmosphere(vec3 p) {
    return exp(-max(p.y, 1) / atmosphereFalloffDepth);
}

vec4 raymarch(vec3 rayOrigin, vec3 rayDirection, vec3 cameraForward, float offset) {
    vec3 lightDirection = normalize(lightPosition);

    float rayDotCam = dot(rayDirection, cameraForward);
    float rayDotCloudPlane = dot(rayDirection, normalize(rayDirection * vec3(1, 0, 1)));
    float distCamCloudPlane = abs(rayOrigin.y - cloudHeight);

    float opaqueDepth = 0;
    vec4 opaqueRes = vec4(0.0);

    for(int i = 0; i < MAX_STEPS; i++) {
        vec3 p = rayOrigin + rayDirection * opaqueDepth;
        float distanceToWorld = sdf(p);

        if(distanceToWorld <= 0.01f) {
            vec3 normal = calculateNormal(p);
            float diffuseIntensity = clamp(dot(lightDirection, normal), 0, 1);
            opaqueRes = mix(vec4(0, 0, 0, 1), vec4(0.2f, 1.0f, 0.6f, 1.0f), diffuseIntensity);
            opaqueRes.rgb *= pointLightColor * pointLightIntensityMultiplier;
            opaqueRes = pow(opaqueRes, vec4(1 / 2.2f)); // Gamma correction
            break;
        }

        opaqueDepth += distanceToWorld;

        if(opaqueDepth > MAX_MARCH_DISTANCE) {
            break;
        }
    }

    vec3 cloudBoxMin = vec3(-cloudBoxWidth, cloudHeight - cloudDepth, -cloudBoxWidth);
    vec3 cloudBoxMax = vec3(cloudBoxWidth, cloudHeight + cloudDepth, cloudBoxWidth);

    float depth = 0;
    vec4 volumetricRes = vec4(0.0);
    float tEnter, tExit;
    if(intersectSphere(rayOrigin, rayDirection, planetOrigin, planetRadius + cloudlessDepth + cloudDepth, tEnter, tExit) && tExit > 0) {
        float startDepth = max(tEnter, 0.0);
        float camAlignedPlane = ceil((startDepth * rayDotCam) / MARCH_SIZE) * MARCH_SIZE;
        depth = camAlignedPlane / rayDotCam;  // World-space depth, but camera-aligned
        startDepth = depth;
    
        float depthTraveledThroughMedium = 0;

        for(int i = 0; i < MAX_STEPS; i++) {

            vec3 p = rayOrigin + rayDirection * depth;

            float density = evaluateDensityAt(p);

            if (density > 0.0) {
                float diffuse = clamp((density - evaluateDensityAt(p + 0.3 * lightDirection)) / 0.3, 0.0, 1.0);
                vec3 lin = ambientColor * ambientIntensity + pointLightIntensityMultiplier * pointLightColor * diffuse;
                vec4 color = vec4(mix(vec3(1.0), vec3(0.0), density), density);
                color.rgb *= lin;
                color.rgb *= color.a;
                volumetricRes += color * (1.0 - volumetricRes.a);
                if (volumetricRes.a >= 0.99) {
                    break;
                }
            }

            depthTraveledThroughMedium += MARCH_SIZE;

            // Samples should be taken at higher frequencies when closer to the camera, i.e. when the depth is low
            float depthFactor = clamp(max(depth - samplingIncreaseDepth, 0) / MAX_MARCH_DISTANCE, 0, 1) * clamp(samplingFalloffDistance / max(distCamCloudPlane, 1), 0, 1);
            depth = startDepth + depthTraveledThroughMedium * (samplingIncreaseFactor * depthFactor  + 1);

            if(depth >= tExit || depthTraveledThroughMedium > MAX_MARCH_DISTANCE || depth > opaqueDepth)
                break;
        }
    }

    vec4 res = mix(opaqueRes, volumetricRes, volumetricRes.a);

    int numAtmosphereSamples = 2;
    for (int i = 0; i < numAtmosphereSamples; i++) {
        vec3 p = rayOrigin + rayDirection * i * ATMOSPHERE_MARCH_SIZE;
        float density = sampleAtmosphere(p);
        vec4 color = vec4(density);
        color.rgb *= ambientColor * ambientIntensity;
        color.rgb *= color.a;
        //res += color * (1.0 - res.a);
    }

    return res;
}


void main() {
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
