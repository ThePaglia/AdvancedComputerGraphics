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
const float ambientIntensity = 1.1;
uniform float cloudTime = 1.0f;

// Raymarching
#define MAX_STEPS 100
const float MARCH_SIZE = 0.08;

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

float scene(vec3 p) {
	float plane = p.y + 1.0;
	float sphere = sdSphere(p + vec3(0, -0.5f, 0), 1.0);
	float f = fbm(p);
	float atmosphere = -1 / (1 + p.y * p.y) * 0.01f;
	float density = min(min(plane, sphere) + f, atmosphere);
	return -density;
}

vec4 raymarch(vec3 rayOrigin, vec3 rayDirection) {
	float depth = 0.0;
	vec3 p = rayOrigin + rayDirection * depth;
	vec4 res = vec4(0.0);
	vec3 lightDirection = normalize(lightPosition);

	for(int i = 0; i < MAX_STEPS; i++) {
		float density = scene(p);

		if(density > 0.0) {
			float diffuse = clamp((scene(p) - scene(p + 0.3 * lightDirection)) / 0.3, 0.0, 1.0);
			vec3 lin = ambientColor * ambientIntensity + pointLightIntensityMultiplier * pointLightColor * diffuse;
			vec4 color = vec4(mix(vec3(1.0, 1.0, 1.0), vec3(0.0, 0.0, 0.0), density), density);
			color.rgb *= lin;
			color.rgb *= color.a;
			res += color * (1.0 - res.a);
		}
		depth += MARCH_SIZE;
		p = rayOrigin + rayDirection * depth;
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
	vec4 res = raymarch(rayOrigin, rayDirection);
	color = res.rgb;

	fragmentColor = vec4(color, 1.0);
}
