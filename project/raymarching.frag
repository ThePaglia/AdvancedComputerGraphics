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
#define MAX_STEPS 1000
#define MAX_MARCH_DISTANCE 1000
const float MARCH_SIZE = 0.05;
uniform sampler2D uBlueNoise;

// Scene parameters
vec3 cloudPosition = vec3(0, 0, -10);
float cloudRadius = 2f;

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
	float f = fbm(p);

    float density = 0;
    for (int i = 0; i < clouds.length; i++) {
        cloud c = clouds[i];
        density = min(density, sdSphere(p + c.position, c.radius) + f);
    }
	
	return -density;
}

float sdf(vec3 p) {
    float dist = 1000;
    for (int i = 0; i < clouds.length; i++) {
        cloud c = clouds[i];
        dist = min(dist, sdSphere(p + c.position, c.radius + 1));
    }
    return dist;
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


vec4 raymarch(vec3 rayOrigin, vec3 rayDirection, vec3 cameraForward, float offset) {
    vec4 res = vec4(0.0);
    float depth = 0.0;
    vec3 lightDirection = normalize(lightPosition);

    float rayDotCam = dot(rayDirection, cameraForward);

    for(int i = 0; i < MAX_STEPS; i++) {
        vec3 p = rayOrigin + rayDirection * depth;
        float dist = sdf(p);

        if(dist < 0.1f) {
            float startDepth = max(depth, 0.0);
            float camPlane = ceil((startDepth * rayDotCam) / MARCH_SIZE) * MARCH_SIZE;
            depth = camPlane / rayDotCam;

            for (int j = i; j < MAX_STEPS; j++) {
                i = j;
                vec3 p = rayOrigin + rayDirection * depth;

                if(sdf(p) > 0.2f)
                    break;

                float density = evaluateDensityAt(p);

                if (density > 0.0) {
                    float diffuse = clamp((density - evaluateDensityAt(p + 0.3 * lightDirection)) / 0.3, 0.0, 1.0);
                    vec3 lin = ambientColor * ambientIntensity + pointLightIntensityMultiplier * pointLightColor * diffuse;
                    vec4 color = vec4(mix(vec3(1.0), vec3(0.0), density), density);
                    color.rgb *= lin;
                    color.rgb *= color.a;
                    res += color * (1.0 - res.a);
                }

                if (res.a >= 0.98) {
                    return res;
                }

                depth += MARCH_SIZE;
            }
        }

        depth += dist;

        if(depth > MAX_MARCH_DISTANCE) {
            break;
        }
    }

    /*
    for(int j = 0; j < clouds.length; j++) {
        cloud c = clouds[j];

        float tEnter, tExit;
        if (!intersectSphere(rayOrigin, rayDirection, c.position, c.radius + 1, tEnter, tExit)) {
            continue;
        }

        float startDepth = max(tEnter, 0.0);
        float camPlane = ceil((startDepth * rayDotCam) / MARCH_SIZE) * MARCH_SIZE;
        float t = camPlane / rayDotCam;

        for (int i = 0; i < MAX_STEPS && t < tExit; i++) {
            vec3 p = rayOrigin + rayDirection * t;
            float density = evaluateDensityAt(p);

            if (density > 0.0) {
                float diffuse = clamp((density - evaluateDensityAt(p + 0.3 * lightDirection)) / 0.3, 0.0, 1.0);
                vec3 lin = ambientColor * ambientIntensity + pointLightIntensityMultiplier * pointLightColor * diffuse;
                vec4 color = vec4(mix(vec3(1.0), vec3(0.0), density), density);
                color.rgb *= lin;
                color.rgb *= color.a;
                res += color * (1.0 - res.a);
            }

            if (res.a >= 0.99) {
                break;
            }

            t += MARCH_SIZE;
        }
    }
    */

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
	vec4 res = raymarch(rayOrigin, rayDirection, uCameraDir, offset);
	color = res.rgb;

	fragmentColor = vec4(color, 1.0);
}
