#version 430

layout (local_size_x = 16, local_size_y = 16) in;

// Output vertex buffer
layout (binding = 0) buffer VertexBuffer {
    vec4 positions[]; // x, y, z, unused
};

layout (binding = 1) buffer ColorBuffer {
    vec4 colors[]; // r, g, b, a
};

// Uniforms from CPU
uniform int latitudes;
uniform int longitudes;
uniform float radius;
uniform float PI = 3.14159265359;

float getHeightOnUnitSphere(vec3 dir) {
    // Example: replace with real 3D noise or fbm
    return 1.0 + 0.2 * sin(10.0 * dir.x) * cos(10.0 * dir.y) * sin(10.0 * dir.z);
}

vec3 getColorForHeight(float h) {
    return vec3(h); // Greyscale for now, or use a better gradient
}

void main() {
    uint i = gl_GlobalInvocationID.y; // latitude index
    uint j = gl_GlobalInvocationID.x; // longitude index

    if (i > uint(latitudes) || j > uint(longitudes)) return;

    float deltaLat = PI / float(latitudes);
    float deltaLon = 2.0 * PI / float(longitudes);
    float latitudeAngle = PI / 2.0 - float(i) * deltaLat;
    float longitudeAngle = float(j) * deltaLon;

    float xy = radius * cos(latitudeAngle);
    float z = radius * sin(latitudeAngle);

    vec3 dir;
    dir.x = xy * cos(longitudeAngle);
    dir.y = xy * sin(longitudeAngle);
    dir.z = z;
    dir = normalize(dir);

    float height = getHeightOnUnitSphere(dir);
    vec3 pos = dir * height;
    vec3 col = getColorForHeight(height);

    uint index = i * uint(longitudes + 1) + j;
    positions[index] = vec4(pos, 1.0);
    colors[index] = vec4(col, 1.0);
}
