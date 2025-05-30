#version 430

layout (local_size_x = 16, local_size_y = 16) in;

// Output vertex buffer
layout (binding = 0) buffer VertexBuffer {
    vec4 positions[]; // x, y, z, unused
};

layout (binding = 1) buffer ColorBuffer {
    vec4 colors[]; // r, g, b, a
};

layout (binding = 2) buffer NormalBuffer {
    vec4 normals[]; // x, y, z, unused
};

// Uniforms from CPU
uniform int latitudes;
uniform int longitudes;
uniform float radius;
uniform float PI = 3.14159265359;
uniform int noiseOctaves = 6;
uniform float noiseLacunarity = 2.0f;
uniform float noiseGain = 0.45f;
uniform float amplitude = 1.0f;
uniform	float frequency = 1.0f;

int Source[] = {
	151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140, 36, 103, 30, 69, 142,
	8, 99, 37, 240, 21, 10, 23, 190, 6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203,
	117, 35, 11, 32, 57, 177, 33, 88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175, 74, 165,
	71, 134, 139, 48, 27, 166, 77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41,
	55, 46, 245, 40, 244, 102, 143, 54, 65, 25, 63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89,
	18, 169, 200, 196, 135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64, 52, 217, 226, 250,
	124, 123, 5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189,
	28, 42, 223, 183, 170, 213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9,
	129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104, 218, 246, 97, 228, 251, 34,
	242, 193, 238, 210, 144, 12, 191, 179, 162, 241, 81, 51, 145, 235, 249, 14, 239, 107, 49, 192, 214, 31,
	181, 199, 106, 157, 184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254, 138, 236, 205, 93, 222, 114,
	67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180 };

const int RandomSize = 256;
const double Sqrt3 = 1.7320508075688772935;
const double Sqrt5 = 2.2360679774997896964;
int _random[RandomSize * 2];

/// Skewing and unskewing factors for 2D, 3D and 4D,
/// some of them pre-multiplied.
const double F2 = 0.5 * (Sqrt3 - 1.0);

const double G2 = (3.0 - Sqrt3) / 6.0;
const double G22 = G2 * 2.0 - 1;

const double F3 = 1.0 / 3.0;
const double G3 = 1.0 / 6.0;

const double F4 = (Sqrt5 - 1.0) / 4.0;
const double G4 = (5.0 - Sqrt5) / 20.0;
const double G42 = G4 * 2.0;
const double G43 = G4 * 3.0;
const double G44 = G4 * 4.0 - 1.0;

/// <summary>
/// Gradient vectors for 3D (pointing to mid points of all edges of a unit
/// cube)
/// </summary>
int Grad3[12][3] =
{
	{1, 1, 0}, {-1, 1, 0}, {1, -1, 0}, {-1, -1, 0}, {1, 0, 1}, {-1, 0, 1}, {1, 0, -1}, {-1, 0, -1}, {0, 1, 1}, {0, -1, 1}, {0, 1, -1}, {0, -1, -1} };

void Randomize(int seed)
{
	if (seed != 0)
	{
		// Shuffle the array using the given seed
		// Unpack the seed into 4 bytes then perform a bitwise XOR operation
		// with each byte
		int shuffleBuffer[4];
		shuffleBuffer[0] = (seed & 0x00ff);
		shuffleBuffer[1] = ((seed & 0xff00) >> 8);
		shuffleBuffer[2] = ((seed & 0x00ff0000) >> 16);
		shuffleBuffer[3] = ((seed & 0xff000000) >> 24);

		for (int i = 0; i < 256; i++)
		{
			_random[i] = Source[i] ^ shuffleBuffer[0];
			_random[i] ^= shuffleBuffer[1];
			_random[i] ^= shuffleBuffer[2];
			_random[i] ^= shuffleBuffer[3];

			_random[i + RandomSize] = _random[i];
		}
	}
	else
	{
		for (int i = 0; i < RandomSize; i++)
			_random[i + RandomSize] = _random[i] = Source[i];
	}
}

float get3DSimplexNoiseAtPointOnUnitSphere(vec3 point)
{
	double x = point.x;
	double y = point.y;
	double z = point.z;
	double n0 = 0, n1 = 0, n2 = 0, n3 = 0;

	// Noise contributions from the four corners
	// Skew the input space to determine which simplex cell we're in
	double s = (x + y + z) * F3;

	// for 3D
	int i = int(floor(x + s));
	int j = int(floor(y + s));
	int k = int(floor(z + s));

	double t = (i + j + k) * G3;

	// The x,y,z distances from the cell origin
	double x0 = x - (i - t);
	double y0 = y - (j - t);
	double z0 = z - (k - t);

	// For the 3D case, the simplex shape is a slightly irregular tetrahedron.
	// Determine which simplex we are in.
	// Offsets for second corner of simplex in (i,j,k)
	int i1, j1, k1;

	// coords
	int i2, j2, k2; // Offsets for third corner of simplex in (i,j,k) coords

	if (x0 >= y0)
	{
		if (y0 >= z0)
		{
			// X Y Z order
			i1 = 1;
			j1 = 0;
			k1 = 0;
			i2 = 1;
			j2 = 1;
			k2 = 0;
		}
		else if (x0 >= z0)
		{
			// X Z Y order
			i1 = 1;
			j1 = 0;
			k1 = 0;
			i2 = 1;
			j2 = 0;
			k2 = 1;
		}
		else
		{
			// Z X Y order
			i1 = 0;
			j1 = 0;
			k1 = 1;
			i2 = 1;
			j2 = 0;
			k2 = 1;
		}
	}
	else
	{
		// x0 < y0
		if (y0 < z0)
		{
			// Z Y X order
			i1 = 0;
			j1 = 0;
			k1 = 1;
			i2 = 0;
			j2 = 1;
			k2 = 1;
		}
		else if (x0 < z0)
		{
			// Y Z X order
			i1 = 0;
			j1 = 1;
			k1 = 0;
			i2 = 0;
			j2 = 1;
			k2 = 1;
		}
		else
		{
			// Y X Z order
			i1 = 0;
			j1 = 1;
			k1 = 0;
			i2 = 1;
			j2 = 1;
			k2 = 0;
		}
	}

	// A step of (1,0,0) in (i,j,k) means a step of (1-c,-c,-c) in (x,y,z),
	// a step of (0,1,0) in (i,j,k) means a step of (-c,1-c,-c) in (x,y,z),
	// and
	// a step of (0,0,1) in (i,j,k) means a step of (-c,-c,1-c) in (x,y,z),
	// where c = 1/6.

	// Offsets for second corner in (x,y,z) coords
	double x1 = x0 - i1 + G3;
	double y1 = y0 - j1 + G3;
	double z1 = z0 - k1 + G3;

	// Offsets for third corner in (x,y,z)
	double x2 = x0 - i2 + F3;
	double y2 = y0 - j2 + F3;
	double z2 = z0 - k2 + F3;

	// Offsets for last corner in (x,y,z)
	double x3 = x0 - 0.5;
	double y3 = y0 - 0.5;
	double z3 = z0 - 0.5;

	// Work out the hashed gradient indices of the four simplex corners
	int ii = i & 0xff;
	int jj = j & 0xff;
	int kk = k & 0xff;

	// Calculate the contribution from the four corners
	double t0 = 0.6 - x0 * x0 - y0 * y0 - z0 * z0;
	if (t0 > 0)
	{
		t0 *= t0;
		int gi0 = _random[ii + _random[jj + _random[kk]]] % 12;
		n0 = t0 * t0 * dot(vec3(Grad3[gi0][0], Grad3[gi0][1], Grad3[gi0][2]), vec3(x0, y0, z0));
	}

	double t1 = 0.6 - x1 * x1 - y1 * y1 - z1 * z1;
	if (t1 > 0)
	{
		t1 *= t1;
		int gi1 = _random[ii + i1 + _random[jj + j1 + _random[kk + k1]]] % 12;
		n1 = t1 * t1 * dot(vec3(Grad3[gi1][0], Grad3[gi1][1], Grad3[gi1][2]), vec3(x1, y1, z1));
	}

	double t2 = 0.6 - x2 * x2 - y2 * y2 - z2 * z2;
	if (t2 > 0)
	{
		t2 *= t2;
		int gi2 = _random[ii + i2 + _random[jj + j2 + _random[kk + k2]]] % 12;
		n2 = t2 * t2 * dot(vec3(Grad3[gi2][0], Grad3[gi2][1], Grad3[gi2][2]), vec3(x2, y2, z2));
	}

	double t3 = 0.6 - x3 * x3 - y3 * y3 - z3 * z3;
	if (t3 > 0)
	{
		t3 *= t3;
		int gi3 = _random[ii + 1 + _random[jj + 1 + _random[kk + 1]]] % 12;
		n3 = t3 * t3 * dot(vec3(Grad3[gi3][0], Grad3[gi3][1], Grad3[gi3][2]), vec3(x3, y3, z3));
	}

	// Add contributions from each corner to get the final noise value.
	// The result is scaled to stay just inside [-1,1]
	return float((n0 + n1 + n2 + n3) * 32);
}

// Helper: Fractal (octave) noise
float fractalSimplexNoise(vec3 point, int octaves = 5, float lacunarity = 2.0f, float gain = 0.5f)
{
	float a = amplitude;
	float f = frequency;
	float sum = 0.0f;
	float maxSum = 0.0f;

	for (int i = 0; i < octaves; ++i)
	{
		sum += a * abs(get3DSimplexNoiseAtPointOnUnitSphere(point * f));
		maxSum += a;
		a *= gain;
		f *= lacunarity;
	}
	return sum / maxSum;
}

float getHeightOnUnitSphere(vec3 point)
{
	// Fractal noise for more detail
	float noise = fractalSimplexNoise(point, noiseOctaves, noiseLacunarity, noiseGain);

	// Sharpen peaks
	noise = pow(noise, 2.5f);

	// Set a sea level threshold
	float seaLevel = 0.0f;

	// Only allow mountains above sea level
	float mountain = max(noise - seaLevel, 0.0f);

	// Scale mountain height (adjust 0.5f for desired max height)
	float height = 1.0f + mountain * 0.8f;

	return height;
}

// Function to get color based on height
vec3 getColorForHeight(float height)
{
	if (height < 1.0f)				   // Below sea level
		return vec3(0.0f, 0.2f, 0.7f); // Deep blue
	else if (height < 1.1f)
		return vec3(0.9f, 0.85f, 0.6f); // Sand/beach
	else if (height < 1.2f)
		return vec3(0.13f, 0.55f, 0.13f); // Green (grass)
	else if (height < 1.3f)
		return vec3(0.5f, 0.36f, 0.33f); // Brown (mountain)
	else
		return vec3(1.0f, 1.0f, 1.0f); // White (snow)
}

vec3 sphericalToCartesian(float lat, float lon) {
    float xy = radius * cos(lat);
    float z = radius * sin(lat);
    return vec3(xy * cos(lon), xy * sin(lon), z);
}

void main() {
    uint i = gl_GlobalInvocationID.y;
    uint j = gl_GlobalInvocationID.x;

    if (i > uint(latitudes) || j > uint(longitudes)) return;

    float deltaLat = PI / float(latitudes);
    float deltaLon = 2.0 * PI / float(longitudes);

    float lat = PI / 2.0 - float(i) * deltaLat;
    float lon = float(j) * deltaLon;

    vec3 dir = normalize(sphericalToCartesian(lat, lon));
    float h = getHeightOnUnitSphere(dir);
    vec3 pos = dir * h;

    // Central difference offsets
    int iUp = int(clamp(int(i) + 1, 0, latitudes));
    int iDown = int(clamp(int(i) - 1, 0, latitudes));
    int jRight = int(j + 1) % (longitudes + 1);
    int jLeft = int(j + longitudes) % (longitudes + 1); // handles wraparound

    float latUp = PI / 2.0 - float(iUp) * deltaLat;
    float latDown = PI / 2.0 - float(iDown) * deltaLat;
    float lonRight = float(jRight) * deltaLon;
    float lonLeft = float(jLeft) * deltaLon;

    vec3 dirLatUp = normalize(sphericalToCartesian(latUp, lon));
    vec3 dirLatDown = normalize(sphericalToCartesian(latDown, lon));
    vec3 dirLonRight = normalize(sphericalToCartesian(lat, lonRight));
    vec3 dirLonLeft = normalize(sphericalToCartesian(lat, lonLeft));

    vec3 posLatUp = dirLatUp * getHeightOnUnitSphere(dirLatUp);
    vec3 posLatDown = dirLatDown * getHeightOnUnitSphere(dirLatDown);
    vec3 posLonRight = dirLonRight * getHeightOnUnitSphere(dirLonRight);
    vec3 posLonLeft = dirLonLeft * getHeightOnUnitSphere(dirLonLeft);

    vec3 dLat = posLatUp - posLatDown;
    vec3 dLon = posLonRight - posLonLeft;

    vec3 normal = normalize(cross(dLat, dLon));

    uint index = i * (longitudes + 1u) + j;
    positions[index] = vec4(pos, h);
    colors[index] = vec4(getColorForHeight(h), 1.0);
    normals[index] = vec4(normal, 1.0);
}