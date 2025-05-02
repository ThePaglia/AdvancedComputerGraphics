#version 420

layout(location = 0) out vec4 fragmentColor;

uniform float uTime;
uniform vec2 uResolution;
// Camera
uniform vec3 uCameraPos;
uniform vec3 uCameraRight;
uniform vec3 uCameraUp;
uniform vec3 uCameraDir;

#define MAX_STEPS 100

float sdSphere(vec3 p, float radius) {
  return length(p) - radius;
}

float scene(vec3 p) {
  float distance = sdSphere(p, 1.0);
  return -distance;
}

const float MARCH_SIZE = 0.08;

vec4 raymarch(vec3 rayOrigin, vec3 rayDirection) {
  float depth = 0.0;
  vec3 p = rayOrigin + rayDirection * depth;
  vec4 res = vec4(0.0);

  for(int i = 0; i < MAX_STEPS; i++) {
    float density = scene(p);

    if(density > 0.0) {
      vec4 color = vec4(mix(vec3(1.0, 1.0, 1.0), vec3(0.0, 0.0, 0.0), density), density);
      color.rgb *= color.a;
      res += color * (1.0 - res.a);
    }
    depth += MARCH_SIZE;
    p = rayOrigin + rayDirection * depth;
  }
  return res;
}

void main() {
  mat3 uCameraMatrix = mat3(uCameraRight, uCameraUp, -uCameraDir);
  vec2 uv = gl_FragCoord.xy / uResolution.xy;
  uv -= 0.5;
  uv.x *= uResolution.x / uResolution.y;

  // Light Position
  vec3 lightPosition = vec3(-10.0, 10.0, 10.0);

  // Ray Origin - camera
  vec3 rayOrigin = uCameraPos;
  // Ray Direction
  vec3 rayDirection = normalize(vec3(uv, -1.0) * uCameraMatrix);
  
  vec3 color = vec3(0.0);
  vec4 res = raymarch(rayOrigin, rayDirection);
  color = res.rgb;
  
  fragmentColor = vec4(color, 1.0);
}
