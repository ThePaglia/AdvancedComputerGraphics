#ifdef _WIN32
extern "C" _declspec(dllexport) unsigned int NvOptimusEnablement = 0x00000001;
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <GL/glew.h>

#include <imgui.h>
#include <labhelper.h>
#include <Model.h>

#include <perf.h>

#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>
using namespace glm;

#include "hdr.h"
#include "fbo.h"
#include <stb_image.h>
#include <iostream>

// Various globals
SDL_Window *g_window = nullptr;
float currentTime = 0.0f;
float previousTime = 0.0f;
float deltaTime = 0.0f;
int windowWidth, windowHeight;

// Mouse input
bool g_doMouseLookaround = false;

///////////////////////////////////////////////////////////////////////////////
// Shader programs
///////////////////////////////////////////////////////////////////////////////
GLuint raymarchingProgram; // Shader for rendering all ray marched effects (clouds, cloud shadows & atmosphere)
GLuint shaderProgram; // Shader for rendering the final image
GLuint depthProgram;  // Shader used to draw the shadow map
GLuint simpleProgram; // Shader used for testing

// Camera parameters.
vec3 worldUp(0.0f, 1.0f, 0.0f);
vec3 cameraPosition(0.0f, 0.0f, -50);
vec3 cameraDirection = normalize(vec3(0.0f) + cameraPosition);
vec3 cameraRight = cross(cameraDirection, worldUp);
vec3 cameraUp = cross(cameraRight, cameraDirection);
mat4 viewProjMatrix;

// Model parameters
labhelper::Model *landingpadModel = nullptr; // Used for debugging the light source's depth buffer
labhelper::Model *sphereModel = nullptr;	 // Used for debug rendering the light source
labhelper::Model *shipModel = nullptr;

float cameraSpeed = 10;

// Texture parameters
GLuint whiteNoiseTexture;

///////////////////////////////////////////////////////////////////////////////
// Light source
///////////////////////////////////////////////////////////////////////////////
vec3 lightPosition = vec3(20, 40, 0);
bool animateLight = false;
vec3 directionalLightColor = vec3(1.0f);
float directionalLightIntensityMultiplier = 1.0f;

// Cloud parameters
float cloudMovementSpeed = 0.05f;
float cloudTime = 0.0f;
float planetRadius = 15.0f;
float cloudlessDepth = 2.2f;
float cloudDepth = 1.0f;
float cloudScale = 0.4f;
float cloudStepMin = 0.01f;
float cloudStepMax = 0.46f;
float cloudShadowIntensity = 2.5f;
float cloudShadowCutoff = 0.4f;
float cloudLightingFalloff = 0.4f;
float cloudNoiseUVScale = 128.0f;
float cloudNoiseAmount = 0.1f;
float sunsetCloudWidth = 0.1f;
int cloudIterations = 6;
int cloudShadowIterations = 4;
float mieIntensity = 0.32f;
float mieG = 0.76f;
bool beerPowderLaw = false;

// Atmosphere parameters
float atmosphereDepth = 10.0f;
float atmosphereDensityFalloff = 3.0f;
vec3 colorBandWavelengths = vec3(700, 530, 440);
float atmosphereScatteringStrength = 3.0f;
float atmosphereDensityAtSeaLevel = 0.17f;
float pointLightIntensityMultiplier = 0.8f;

// Shadow map settings
FboInfo shadowMapFB;
int shadowMapResolution = 4096;
bool usePolygonOffset = true;
bool useHardwarePCF = true;
float polygonOffset_factor = 2.0f;
float polygonOffset_units = 30000.0f;

// Render textures, these are rendered to in the rasterizer step, and the result is used in the ray marching shader
FboInfo rasterizedFBO;
GLuint colorTex;
GLuint depthTex;

// Planet buffers etc.
GLuint planetPositionBuffer, planetIndexBuffer, planetNormalBuffer, planetVAO;
std::vector<vec3> planetVertices;
std::vector<vec3> planetNormals;
std::vector<vec2> planetUVs;
std::vector<unsigned int> planetIndices;
mat4 planetModelMatrix;



float getHeightOnUnitSphere(vec3 p) 
{
	return 1.0f - sin(p.x * 10) * 0.1f;
}

// Modified code from https://gist.github.com/Pikachuxxxx/5c4c490a7d7679824e0e18af42918efc
void generateSphereSmooth(int radius, int latitudes, int longitudes)
{
	if (longitudes < 3)
		longitudes = 3;
	if (latitudes < 2)
		latitudes = 2;

	std::vector<vec3> vertices;
	std::vector<vec3> normals;
	std::vector<vec2> uv;
	std::vector<unsigned int> indices;

	float nx, ny, nz, lengthInv = 1.0f / radius;    // normal
	// Temporary vertex
	struct Vertex
	{
		float x, y, z, s, t; // Postion and Texcoords
	};

	float deltaLatitude = M_PI / latitudes;
	float deltaLongitude = 2 * M_PI / longitudes;
	float latitudeAngle;
	float longitudeAngle;

	// Compute all vertices first except normals
	for (int i = 0; i <= latitudes; ++i)
	{
		latitudeAngle = M_PI / 2 - i * deltaLatitude; /* Starting -pi/2 to pi/2 */
		float xy = radius * cosf(latitudeAngle);    /* r * cos(phi) */
		float z = radius * sinf(latitudeAngle);     /* r * sin(phi )*/

		/*
		 * We add (latitudes + 1) vertices per longitude because of equator,
		 * the North pole and South pole are not counted here, as they overlap.
		 * The first and last vertices have same position and normal, but
		 * different tex coords.
		 */
		for (int j = 0; j <= longitudes; ++j)
		{
			longitudeAngle = j * deltaLongitude;

			// dir points from the center of the sphere towards this vertex's point on the unit sphere. It is used to scale the points according to noise
			vec3 dir;
			dir.x = xy * cosf(longitudeAngle);       /* x = r * cos(phi) * cos(theta)  */
			dir.y = xy * sinf(longitudeAngle);       /* y = r * cos(phi) * sin(theta) */
			dir.z = z;                               /* z = r * sin(phi) */
			dir = normalize(dir);

			// The terrain height in this direction
			float height = getHeightOnUnitSphere(dir);
			
			// The final vertex position is acquired by multiplying the dir by the height value, we effectively get a heightmap mapped onto the sphere
			Vertex vertex;
			vertex.x = dir.x * height;
			vertex.y = dir.y * height;
			vertex.z = dir.z * height;
			vertex.s = (float)j / longitudes;             /* s */
			vertex.t = (float)i / latitudes;              /* t */
			vertices.push_back(glm::vec3(vertex.x, vertex.y, vertex.z));
			uv.push_back(glm::vec2(vertex.s, vertex.t));

			// normalized vertex normal
			nx = vertex.x * lengthInv;
			ny = vertex.y * lengthInv;
			nz = vertex.z * lengthInv;
			normals.push_back(glm::vec3(nx, ny, nz));
		}
	}

	/*
	 *  Indices
	 *  k1--k1+1
	 *  |  / |
	 *  | /  |
	 *  k2--k2+1
	 */
	unsigned int k1, k2;
	for (int i = 0; i < latitudes; ++i)
	{
		k1 = i * (longitudes + 1);
		k2 = k1 + longitudes + 1;
		// 2 Triangles per latitude block excluding the first and last longitudes blocks
		for (int j = 0; j < longitudes; ++j, ++k1, ++k2)
		{
			if (i != 0)
			{
				indices.push_back(k1);
				indices.push_back(k2);
				indices.push_back(k1 + 1);
			}

			if (i != (latitudes - 1))
			{
				indices.push_back(k1 + 1);
				indices.push_back(k2);
				indices.push_back(k2 + 1);
			}
		}
	}

	planetVertices = vertices;
	planetNormals = normals;
	planetUVs = uv;
	planetIndices = indices;
}

void initializePlanet()
{
	generateSphereSmooth(1, 100, 100);

	///////////////////////////////////////////////////////////////////////////
	// Create the vertex array object
	///////////////////////////////////////////////////////////////////////////
	// Create a handle for the vertex array object
	glGenVertexArrays(1, &planetVAO);
	// Set it as current, i.e., related calls will affect this object
	glBindVertexArray(planetVAO);

	// Create a handle for the vertex position buffer
	glGenBuffers(1, &planetPositionBuffer);
	// Set the newly created buffer as the current one
	glBindBuffer(GL_ARRAY_BUFFER, planetPositionBuffer);
	// Send the vertex position data to the current buffer
	glBufferData(GL_ARRAY_BUFFER, planetVertices.size() * sizeof(vec3), planetVertices.data(), GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, false /*normalized*/, 0 /*stride*/, 0 /*offset*/);
	// Enable the attribute
	glEnableVertexAttribArray(0);

	// Create a handle for the normal buffer
	glGenBuffers(1, &planetNormalBuffer);
	// Set the newly created buffer as the current one
	glBindBuffer(GL_ARRAY_BUFFER, planetNormalBuffer);
	// Send the vertex normal data to the current buffer
	glBufferData(GL_ARRAY_BUFFER, planetNormals.size() * sizeof(vec3), planetNormals.data(), GL_STATIC_DRAW);
	glVertexAttribPointer(1, 3, GL_FLOAT, false /*normalized*/, 0 /*stride*/, 0 /*offset*/);
	// Enable the attribute
	glEnableVertexAttribArray(1);

	// Finally, set the planet's indices
	glGenBuffers(1, &planetIndexBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, planetIndexBuffer);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, planetIndices.size() * sizeof(unsigned int), planetIndices.data(),
		GL_STATIC_DRAW);
}

void loadModels() {
	landingpadModel = labhelper::loadModelFromOBJ("../scenes/landingpad.obj");
	sphereModel = labhelper::loadModelFromOBJ("../scenes/sphere.obj");
	shipModel = labhelper::loadModelFromOBJ("../scenes/NewShip.obj");
}

void loadShaders(bool is_reload)
{
	GLuint shader = labhelper::loadShaderProgram("../project/depth.vert", "../project/depth.frag", is_reload);
	if (shader != 0)
	{
		depthProgram = shader;
	}

	shader = labhelper::loadShaderProgram("../project/shading.vert", "../project/shading.frag", is_reload);
	if (shader != 0)
	{
		shaderProgram = shader;
	}

	shader = labhelper::loadShaderProgram("../project/raymarching.vert", "../project/raymarching.frag", is_reload);
	if (shader != 0)
	{
		raymarchingProgram = shader;
	}

	shader = labhelper::loadShaderProgram("../project/simple.vert", "../project/simple.frag", is_reload);
	if (shader != 0)
	{
		simpleProgram = shader;
	}
}

void loadNoiseTexture(const std::string& filepath, GLuint& texture)
{
	int width, height, channels;
	unsigned char* data = stbi_load(filepath.c_str(), &width, &height, &channels, 0);
	if (!data)
	{
		std::cerr << "Failed to load texture: " << filepath << std::endl;
		return;
	}

	glGenTextures(1, &texture);
	glBindTexture(GL_TEXTURE_2D, texture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
	glGenerateMipmap(GL_TEXTURE_2D);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	stbi_image_free(data);
}


// This function is called once at the start of the program and never again
void initialize()
{
	ENSURE_INITIALIZE_ONLY_ONCE();

	// Load Shaders
	loadShaders(false);

	// Load models
	loadModels();

	// Initialize the planet
	initializePlanet();

	// Load noise texture
	loadNoiseTexture("../textures/whiteNoise.png", whiteNoiseTexture);

	// Shadow map
	shadowMapFB.resize(shadowMapResolution, shadowMapResolution);
	glBindTexture(GL_TEXTURE_2D, shadowMapFB.depthBuffer);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_COMPARE_REF_TO_TEXTURE);

	glEnable(GL_DEPTH_TEST); // enable Z-buffering
	glEnable(GL_CULL_FACE);	 // enables backface culling
}

// This function is used to draw the main objects on the scene
void drawFullscreenQuad(GLuint currentShaderProgram)
{
	glUseProgram(currentShaderProgram);
	glFrontFace(GL_CCW); // The drawing order is flipped for the full-screen quad used for raymarching

	// Bind the noise texture
	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_2D, whiteNoiseTexture);
	labhelper::setUniformSlow(currentShaderProgram, "uWhiteNoiseTexture", 2);

	// Light source
	labhelper::setUniformSlow(currentShaderProgram, "directionalLightColor", directionalLightColor);
	labhelper::setUniformSlow(currentShaderProgram, "directionalLightIntensityMultiplier",
							  directionalLightIntensityMultiplier);
	labhelper::setUniformSlow(currentShaderProgram, "lightPosition", lightPosition);

	// Simulation parameters
	labhelper::setUniformSlow(currentShaderProgram, "uTime", currentTime);
	labhelper::setUniformSlow(currentShaderProgram, "cloudTime", cloudTime);
	labhelper::setUniformSlow(currentShaderProgram, "planetRadius", planetRadius);
	labhelper::setUniformSlow(currentShaderProgram, "cloudlessDepth", cloudlessDepth);
	labhelper::setUniformSlow(currentShaderProgram, "cloudDepth", cloudDepth);
	labhelper::setUniformSlow(currentShaderProgram, "cloudScale", cloudScale);
	labhelper::setUniformSlow(currentShaderProgram, "cloudStepMin", cloudStepMin);
	labhelper::setUniformSlow(currentShaderProgram, "cloudShadowCutoff", cloudShadowCutoff);
	labhelper::setUniformSlow(currentShaderProgram, "cloudShadowIntensity", cloudShadowIntensity);
	labhelper::setUniformSlow(currentShaderProgram, "cloudStepMax", cloudStepMax);
	labhelper::setUniformSlow(currentShaderProgram, "cloudLightingFalloff", cloudLightingFalloff);
	labhelper::setUniformSlow(currentShaderProgram, "cloudNoiseUVScale", cloudNoiseUVScale);
	labhelper::setUniformSlow(currentShaderProgram, "cloudNoiseAmount", cloudNoiseAmount);
	labhelper::setUniformSlow(currentShaderProgram, "cloudShadowIterations", cloudShadowIterations);
	labhelper::setUniformSlow(currentShaderProgram, "sunsetCloudWidth", sunsetCloudWidth);
	labhelper::setUniformSlow(currentShaderProgram, "cloudIterations", cloudIterations);
	labhelper::setUniformSlow(currentShaderProgram, "atmosphereDepth", atmosphereDepth);
	labhelper::setUniformSlow(currentShaderProgram, "atmosphereDensityFalloff", atmosphereDensityFalloff);
	labhelper::setUniformSlow(currentShaderProgram, "atmosphereDensityAtSeaLevel", atmosphereDensityAtSeaLevel);
	vec3 scatteringCoefficients = vec3(pow(300 / colorBandWavelengths.x, 4), pow(300 / colorBandWavelengths.y, 4), pow(300 / colorBandWavelengths.z, 4)) * atmosphereScatteringStrength;
	labhelper::setUniformSlow(currentShaderProgram, "atmosphereScatteringCoefficients", scatteringCoefficients);
	labhelper::setUniformSlow(currentShaderProgram, "mieIntensity", mieIntensity);
	labhelper::setUniformSlow(currentShaderProgram, "mieG", mieG);
	labhelper::setUniformSlow(currentShaderProgram, "beerPowderLaw", beerPowderLaw);

	// uResolution
	labhelper::setUniformSlow(currentShaderProgram, "uResolution", vec2(windowWidth, windowHeight));

	// Camera
	labhelper::setUniformSlow(currentShaderProgram, "uCameraPos", cameraPosition);
	labhelper::setUniformSlow(currentShaderProgram, "uCameraDir", cameraDirection);
	labhelper::setUniformSlow(currentShaderProgram, "uCameraUp", cameraUp);
	labhelper::setUniformSlow(currentShaderProgram, "uCameraRight", cameraRight);
	labhelper::setUniformSlow(currentShaderProgram, "uViewProjectionMatrix", viewProjMatrix);

	labhelper::drawFullScreenQuad();
}

void drawScene(GLuint currentShaderProgram,
					   const mat4 &viewMatrix,
					   const mat4 &projectionMatrix,
					   const mat4 &lightViewMatrix,
					   const mat4 &lightProjectionMatrix)
{
	glUseProgram(currentShaderProgram);
	if (currentShaderProgram == depthProgram)
	{
		glFrontFace(GL_CCW); // depthProgram requires CCW vertex order for some reason, this way it properly renders the forward facing faces
	}
	else
	{
		glFrontFace(GL_CW); // The models are rendered inside out so we flip what is considered to be the front face
	}

	// Light source
	vec4 viewSpaceLightPosition = viewMatrix * vec4(lightPosition, 1.0f);
	labhelper::setUniformSlow(currentShaderProgram, "directional_light_color", directionalLightColor);
	labhelper::setUniformSlow(currentShaderProgram, "directional_light_intensity_multiplier", directionalLightIntensityMultiplier);
	labhelper::setUniformSlow(currentShaderProgram, "viewSpaceLightPosition", vec3(viewSpaceLightPosition));
	labhelper::setUniformSlow(currentShaderProgram, "viewSpaceLightDir", normalize(vec3(viewMatrix * vec4(-lightPosition, 0.0f))));

	// Camera
	labhelper::setUniformSlow(currentShaderProgram, "viewInverse", inverse(viewMatrix));

	// Light matrix and depth buffer texture
	glActiveTexture(GL_TEXTURE10);
	glBindTexture(GL_TEXTURE_2D, shadowMapFB.depthBuffer);
	mat4 lightMatrix = translate(vec3(0.5f)) * scale(vec3(0.5f)) * lightProjectionMatrix * lightViewMatrix * inverse(viewMatrix);
	labhelper::setUniformSlow(currentShaderProgram, "lightMatrix", lightMatrix);

	// Uncomment this to render the landing pad, used for debugging the light's depthBuffer
	/*
	mat4 modelMatrix = translate(vec3(0, -30, 0)) * rotate(radians(-90.0f), worldUp);
	labhelper::setUniformSlow(currentShaderProgram, "modelViewProjectionMatrix",
		projectionMatrix * viewMatrix * modelMatrix);
	labhelper::setUniformSlow(currentShaderProgram, "modelViewMatrix", viewMatrix * modelMatrix);
	labhelper::setUniformSlow(currentShaderProgram, "normalMatrix",
		inverse(transpose(viewMatrix * modelMatrix)));
	labhelper::render(landingpadModel);
	*/


	// Planet
	labhelper::setUniformSlow(currentShaderProgram, "modelViewProjectionMatrix", projectionMatrix * viewMatrix * planetModelMatrix);
	labhelper::setUniformSlow(currentShaderProgram, "modelViewMatrix", viewMatrix * planetModelMatrix);
	labhelper::setUniformSlow(currentShaderProgram, "normalMatrix", inverse(transpose(viewMatrix * planetModelMatrix)));
	labhelper::setUniformSlow(currentShaderProgram, "material_color", vec3(0.1f, 1.0f, 0.3f));
	glBindVertexArray(planetVAO);
	glDrawElements(GL_TRIANGLES, planetIndices.size(), GL_UNSIGNED_INT, 0);

	// Ship
	float d = 16 + sin(cloudTime * 2) * 0.25f;
	mat4 shipMatrix = rotate(radians(25.0f), vec3(0, 0, 1)) * translate(vec3(0, d, 0)) * scale(vec3(0.05f));
	labhelper::setUniformSlow(currentShaderProgram, "modelViewProjectionMatrix", projectionMatrix * viewMatrix * shipMatrix);
	labhelper::setUniformSlow(currentShaderProgram, "modelViewMatrix", viewMatrix * shipMatrix);
	labhelper::setUniformSlow(currentShaderProgram, "normalMatrix", inverse(transpose(viewMatrix * shipMatrix)));
	labhelper::render(shipModel);
}

// This function will be called once per frame, so the code to set up
// the scene for rendering should go here
void display(void)
{
	labhelper::perf::Scope s("Display");

	// Check if window size has changed and resize buffers as needed
	{
		int w, h;
		SDL_GetWindowSize(g_window, &w, &h);
		if (w != windowWidth || h != windowHeight)
		{
			windowWidth = w;
			windowHeight = h;
			rasterizedFBO.resize(w, h);
		}
	}

	// setup matrices
	mat4 projMatrix = perspective(radians(45.0f), float(windowWidth) / float(windowHeight), 0.01f, 1000.0f);
	mat3 cameraBaseVectorsWorldSpace(cameraRight, cameraUp, cameraDirection);
	mat4 cameraRotation = mat4(transpose(cameraBaseVectorsWorldSpace)); // NOTE: this is also calculated in the raymarching shader, perhaps we can just send the result there?
	mat4 viewMatrix = cameraRotation * translate(-cameraPosition);
	viewProjMatrix = projMatrix * viewMatrix;
	// viewMatrix = lookAt(cameraPosition, cameraPosition - cameraDirection, worldUp);

	vec4 lightStartPosition = vec4(40.0f, 40.0f, 0.0f, 1.0f);
	if (animateLight)
		lightPosition = vec3(rotate(currentTime * 0.5f, worldUp) * lightStartPosition);

	mat4 lightViewMatrix = lookAt(lightPosition, vec3(0.0f), worldUp);
	// We scale the orthographic "frustum" to the planet's radius so that we get the most out of the shadow map
	mat4 lightProjMatrix = ortho(-planetRadius, planetRadius, -planetRadius, planetRadius, 10.0f, 100.0f);

	///////////////////////////////////////////////////////////////////////////
	// Set Up and draw to Shadow Map
	///////////////////////////////////////////////////////////////////////////
	{
		labhelper::perf::Scope shadowMapScope("Shadow Map");

		if (shadowMapFB.width != shadowMapResolution || shadowMapFB.height != shadowMapResolution)
		{
			shadowMapFB.resize(shadowMapResolution, shadowMapResolution);
		}

		glBindTexture(GL_TEXTURE_2D, shadowMapFB.depthBuffer);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

		if (useHardwarePCF)
		{
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		}
		else
		{
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		}

		// This line is to avoid some warnings from OpenGL for having the shadowmap attached to texture unit 0
		// when using a shader that samples from that texture with a sampler2D instead of a shadow sampler.
		// It is never actually sampled, but just having it set there generates the warning in some systems.
		glBindTexture(GL_TEXTURE_2D, 0);

		///////////////////////////////////////////////////////////////////////////
		// Draw Shadow Map
		///////////////////////////////////////////////////////////////////////////
		if (usePolygonOffset)
		{
			glEnable(GL_POLYGON_OFFSET_FILL);
			glPolygonOffset(polygonOffset_factor, polygonOffset_units);
		}

		glBindTexture(GL_TEXTURE_2D, shadowMapFB.depthBuffer);
		glBindFramebuffer(GL_FRAMEBUFFER, shadowMapFB.framebufferId);
		glViewport(0, 0, shadowMapFB.width, shadowMapFB.height);
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glClearDepth(1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		drawScene(depthProgram, lightViewMatrix, lightProjMatrix, lightViewMatrix, lightProjMatrix);
		if (usePolygonOffset)
		{
			glDisable(GL_POLYGON_OFFSET_FILL);
		}
	}

	// These following two lines are used for debugging the ligth's depthBuffer
	// labhelper::Material& screen = landingpadModel->m_materials[8];
	// screen.m_emission_texture.gl_id = shadowMapFB.colorTextureTargets[0];

	///////////////////////////////////////////////////////////////////////////
	// Set Up and draw to Rasterized FBO texture
	///////////////////////////////////////////////////////////////////////////
	{
		labhelper::perf::Scope rasterizationScope("Rasterized Graphics");

		// Draw to fbo from camera
		glBindFramebuffer(GL_FRAMEBUFFER, rasterizedFBO.framebufferId);
		glViewport(0, 0, rasterizedFBO.width, rasterizedFBO.height);
		glClearColor(0, 0, 0, 0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		planetModelMatrix = scale(vec3(planetRadius));
		drawScene(shaderProgram, viewMatrix, projMatrix, lightViewMatrix, lightProjMatrix);
	}

	///////////////////////////////////////////////////////////////////////////
	// Set Up and draw post process Ray Marching 
	///////////////////////////////////////////////////////////////////////////
	{
		labhelper::perf::Scope rayMarchingScope("Ray Marching");
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// Bind the raymarch shader
		glUseProgram(raymarchingProgram);
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, rasterizedFBO.colorTextureTargets[0]);
		glUniform1i(glGetUniformLocation(raymarchingProgram, "uSceneColor"), 0);

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, rasterizedFBO.depthBuffer);
		glUniform1i(glGetUniformLocation(raymarchingProgram, "uSceneDepth"), 1);

		drawFullscreenQuad(raymarchingProgram);
	}
}

// This function is used to update the scene according to user input
bool handleEvents(void)
{
	// check events (keyboard among other)
	SDL_Event event;
	bool quitEvent = false;
	while (SDL_PollEvent(&event))
	{
		labhelper::processEvent(&event);

		if (event.type == SDL_QUIT || (event.type == SDL_KEYUP && event.key.keysym.sym == SDLK_ESCAPE))
		{
			quitEvent = true;
		}
		if (event.type == SDL_KEYUP && event.key.keysym.sym == SDLK_g)
		{
			if (labhelper::isGUIvisible())
			{
				labhelper::hideGUI();
			}
			else
			{
				labhelper::showGUI();
			}
		}
		if (event.type == SDL_MOUSEBUTTONDOWN && event.button.button == SDL_BUTTON_LEFT && (!labhelper::isGUIvisible() || !ImGui::GetIO().WantCaptureMouse))
		{
			int x;
			int y;
			SDL_GetMouseState(&x, &y);
		}
		if (event.type == SDL_MOUSEBUTTONUP && event.button.button == SDL_BUTTON_LEFT && (!labhelper::isGUIvisible() || !ImGui::GetIO().WantCaptureMouse))
		{
			g_doMouseLookaround = !g_doMouseLookaround;
			SDL_SetRelativeMouseMode(g_doMouseLookaround ? SDL_TRUE : SDL_FALSE);
		}

		if (event.type == SDL_MOUSEMOTION && g_doMouseLookaround)
		{
			// More info at https://wiki.libsdl.org/SDL_MouseMotionEvent
			int delta_x = event.motion.xrel;
			int delta_y = event.motion.yrel;
			float rotationSpeed = -0.1f;
			// Calculate yaw, i.e. rotation around y axis
			mat4 yaw = rotate(rotationSpeed * deltaTime * -delta_x, worldUp);
			// Apply yaw to direction
			cameraDirection = vec3(yaw * vec4(cameraDirection, 0.0f));
			// Calculate right vector from new direction
			cameraRight = normalize(cross(cameraDirection, worldUp));
			// Calculate pitch around new cameraRight
			mat4 pitch = rotate(rotationSpeed * deltaTime * -delta_y, cameraRight);
			// Apply pitch to direction
			cameraDirection = vec3(pitch * vec4(cameraDirection, 0.0f));
			// Calculate cameraUp from new direction and cameraRight
			cameraUp = cross(cameraRight, cameraDirection);
		}
	}

	// check keyboard state (which keys are still pressed)
	const uint8_t *state = SDL_GetKeyboardState(nullptr);

	if (state[SDL_SCANCODE_W])
	{
		cameraPosition -= cameraSpeed * deltaTime * cameraDirection;
	}
	if (state[SDL_SCANCODE_S])
	{
		cameraPosition += cameraSpeed * deltaTime * cameraDirection;
	}
	if (state[SDL_SCANCODE_A])
	{
		cameraPosition -= cameraSpeed * deltaTime * cameraRight;
	}
	if (state[SDL_SCANCODE_D])
	{
		cameraPosition += cameraSpeed * deltaTime * cameraRight;
	}
	if (state[SDL_SCANCODE_LCTRL])
	{
		cameraPosition -= cameraSpeed * deltaTime * worldUp;
	}
	if (state[SDL_SCANCODE_LSHIFT])
	{
		cameraPosition += cameraSpeed * deltaTime * worldUp;
	}
	if (state[SDL_SCANCODE_E])
	{
		// cameraPosition += cameraSpeed * deltaTime * normalize(cross(cameraDirection, worldUp));
	}
	if (state[SDL_SCANCODE_Q])
	{
		// cameraPosition -= cameraSpeed * deltaTime * normalize(cross(cameraDirection, worldUp));
	}

	return quitEvent;
}

// This function is to hold the general GUI logic

void gui()
{
	// ----------------- Set variables --------------------------
	ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate,
				ImGui::GetIO().Framerate);
	// ----------------------------------------------------------
	ImGui::Text("Settings");
	ImGui::SliderFloat("Camera Speed", &cameraSpeed, 0.1f, 100.0f);
	ImGui::SliderFloat("Planet Radius", &planetRadius, 1, 100.0f);

	ImGui::Text("Clouds");
	ImGui::SliderInt("Cloud Iterations", &cloudIterations, 1, 10);
	ImGui::SliderInt("Cloud Shadow Iterations", &cloudShadowIterations, 1, 6);
	ImGui::SliderFloat("Cloud Speed", &cloudMovementSpeed, 0.0f, 1.0f);
	ImGui::SliderFloat("Cloudless Depth", &cloudlessDepth, 0.1f, 10);
	ImGui::SliderFloat("Cloud Depth", &cloudDepth, 0.0f, 1);
	ImGui::SliderFloat("Cloud Scale", &cloudScale, 0.1f, 1);
	ImGui::SliderFloat("Cloud Step Min", &cloudStepMin, 0, 1);
	ImGui::SliderFloat("Cloud Step Max", &cloudStepMax, 0, 1);
	ImGui::SliderFloat("Cloud Shadow Cutoff", &cloudShadowCutoff, 0, 1);
	ImGui::SliderFloat("Cloud Shadow Intensity", &cloudShadowIntensity, 0, 10);
	ImGui::SliderFloat("Cloud Lighting Fraction", &cloudLightingFalloff, 0.001f, 1.0f);
	ImGui::SliderFloat("Cloud Noise UV Scale", &cloudNoiseUVScale, 1, 2048);
	ImGui::SliderFloat("Cloud Noise Amount", &cloudNoiseAmount, 0, 1);
	ImGui::SliderFloat("Sunset Cloud Width", &sunsetCloudWidth, 0.001f, 1.0f);
	ImGui::SliderFloat("Mie Scattering Intensity", &mieIntensity, 0.0, 1.0);
	ImGui::SliderFloat("Mie Scattering G", &mieG, 0.0f, 0.999f);
	ImGui::Checkbox("Beer-Powder Law", &beerPowderLaw);
	ImGui::Text("Atmosphere");
	ImGui::SliderFloat("Atmosphere Depth", &atmosphereDepth, 0, 10);
	ImGui::SliderFloat("Density Falloff", &atmosphereDensityFalloff, 0, 10);
	ImGui::SliderFloat("Density at Sea Level", &atmosphereDensityAtSeaLevel, 0, 1);
	ImGui::SliderFloat3("Scattering Wavelengths", (float *)&colorBandWavelengths, 0, 1000);
	ImGui::SliderFloat("Scattering Strength", &atmosphereScatteringStrength, 0, 10);

	ImGui::Text("Sun");
	ImGui::SliderFloat("Sun Intensity", &directionalLightIntensityMultiplier, 0, 2);
	ImGui::Checkbox("Animate light", &animateLight);
	ImGui::SliderFloat3("Sun Position", (float *)&lightPosition, -100, 100);
	ImGui::SliderInt("Shadow Map Resolution", &shadowMapResolution, 32, 4096);
	ImGui::Checkbox("Use polygon offset", &usePolygonOffset);
	ImGui::SliderFloat("Factor", &polygonOffset_factor, 0.0f, 10.0f);
	ImGui::SliderFloat("Units", &polygonOffset_units, 0.0f, 100000.0f);
	ImGui::Checkbox("Use hardware PCF", &useHardwarePCF);

	labhelper::perf::drawEventsWindow();
}

int main(int argc, char *argv[])
{
	g_window = labhelper::init_window_SDL("OpenGL Project");

	initialize();

	bool stopRendering = false;
	auto startTime = std::chrono::system_clock::now();

	while (!stopRendering)
	{
		// update currentTime
		std::chrono::duration<float> timeSinceStart = std::chrono::system_clock::now() - startTime;
		previousTime = currentTime;
		currentTime = timeSinceStart.count();
		deltaTime = currentTime - previousTime;
		cloudTime += deltaTime * cloudMovementSpeed;

		// check events (keyboard among other)
		stopRendering = handleEvents();

		// Inform imgui of new frame
		labhelper::newFrame(g_window);

		// render to window
		display();

		// Render overlay GUI.
		gui();

		// Finish the frame and render the GUI
		labhelper::finishFrame();

		// Swap front and back buffer. This frame will now been displayed.
		SDL_GL_SwapWindow(g_window);
	}

	// Shut down everything. This includes the window and all other subsystems.
	labhelper::shutDown(g_window);
	return 0;
}