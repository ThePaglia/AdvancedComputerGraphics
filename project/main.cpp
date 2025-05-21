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

#include <Model.h>
#include "hdr.h"
#include "fbo.h"
#include <stb_image.h>
#include <iostream>

// Various globals
SDL_Window* g_window = nullptr;
float currentTime = 0.0f;
float previousTime = 0.0f;
float deltaTime = 0.0f;
int windowWidth, windowHeight;

// Mouse input
bool g_doMouseLookaround = false;

// This will identify our VBO for the planet
GLuint vertexArrayObject;

// Shader programs
GLuint raymarchingProgram;

///////////////////////////////////////////////////////////////////////////////
// Shader programs
///////////////////////////////////////////////////////////////////////////////
GLuint shaderProgram;       // Shader for rendering the final image
GLuint depthProgram; // Shader used to draw the shadow map

// Camera parameters.
vec3 worldUp(0.0f, 1.0f, 0.0f);
vec3 cameraPosition(0.0f, 0.0f, -30);
vec3 cameraDirection = normalize(vec3(0.0f) + cameraPosition);
vec3 cameraRight = cross(cameraDirection, worldUp);
vec3 cameraUp = cross(cameraRight, cameraDirection);
mat4 viewProjMatrix;

// Model parameters
labhelper::Model* planetModel = nullptr;
mat4 planetModelMatrix;

float cameraSpeed = 10;

// Texture parameters
GLuint noiseTexture;

///////////////////////////////////////////////////////////////////////////////
// Light source
///////////////////////////////////////////////////////////////////////////////
vec3 lightPosition = vec3(0.0f, 100.0f, 0.0f);
bool animateLight = false;
vec3 pointLightColor = vec3(1.0f);
float innerSpotlightAngle = 17.5f;
float outerSpotlightAngle = 22.5f;
float point_light_intensity_multiplier = 5000.0f;

// Sampling parameters
float samplingIncreaseFactor = 20.0f;
float samplingIncreaseDepth = 5.0f;
float samplingFalloffDistance = 10.0f;

// Scene parameters
float cloudMovementSpeed = 0.1f;
float cloudTime = 0.0f;
float planetRadius = 10.0f;
float cloudlessDepth = 0.5f;
float cloudDepth = 1.0f;
float cloudScale = 0.72f;
float cloudStepMin = 0.01f;
float cloudStepMax = 0.46f;
float cloudShadowIntensity = 3.0f;
float cloudShadowCutoff = 0.1f;
float cloudLightingFalloff = 0.5f;
float atmosphereDepth = 5.6f;
float atmosphereDensityFalloff = 5.5f;
vec3 colorBandWavelengths = vec3(700, 530, 440);
float atmosphereScatteringStrength = 3.0f;
float atmosphereDensityAtSeaLevel = 0.5f;
float pointLightIntensityMultiplier = 0.8f;

// Shadow map
enum ClampMode
{
	Edge = 1,
	Border = 2
};

FboInfo shadowMapFB;
int shadowMapResolution = 1024;
int shadowMapClampMode = ClampMode::Edge;
bool shadowMapClampBorderShadowed = false;
bool usePolygonOffset = true;
bool useSoftFalloff = false;
bool useHardwarePCF = false;
float polygonOffset_factor = 1.0f;
float polygonOffset_units = 5000.0f;

// Render textures, these are rendered to in the rasterizer step, and the result is used in the ray marching shader
FboInfo rasterizedFBO;
GLuint colorTex;
GLuint depthTex;

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
}

void loadNoiseTexture(const std::string& filepath)
{
	int width, height, channels;
	unsigned char* data = stbi_load(filepath.c_str(), &width, &height, &channels, 0);
	if (!data)
	{
		std::cerr << "Failed to load texture: " << filepath << std::endl;
		return;
	}

	glGenTextures(1, &noiseTexture);
	glBindTexture(GL_TEXTURE_2D, noiseTexture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
	glGenerateMipmap(GL_TEXTURE_2D);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	stbi_image_free(data);
}

void initializePlanet() {
	planetModel = labhelper::loadModelFromOBJ("../scenes/planet.obj");
	planetModelMatrix = scale(vec3(9));
}

// This function is called once at the start of the program and never again
void initialize()
{
	ENSURE_INITIALIZE_ONLY_ONCE();
	
	// Load Shaders
	loadShaders(false);

	initializePlanet();

	// Load noise texture
	loadNoiseTexture("../textures/noise.png");

	// Shadow map
	shadowMapFB.resize(shadowMapResolution, shadowMapResolution);
	glBindTexture(GL_TEXTURE_2D, shadowMapFB.depthBuffer);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_COMPARE_REF_TO_TEXTURE);

	glEnable(GL_DEPTH_TEST); // enable Z-buffering
	glEnable(GL_CULL_FACE);  // enables backface culling
}

// This function is used to draw the main objects on the scene
void drawScene(GLuint currentShaderProgram,
	const mat4& viewMatrix,
	const mat4& projectionMatrix)
{
	glUseProgram(currentShaderProgram);
	glFrontFace(GL_CCW); // The drawing order is flipped for the full-screen quad used for raymarching

	// Bind the noise texture
	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_2D, noiseTexture);
	labhelper::setUniformSlow(currentShaderProgram, "uNoiseTexture", 2);

	// Light source
	labhelper::setUniformSlow(currentShaderProgram, "pointLightColor", pointLightColor);
	labhelper::setUniformSlow(currentShaderProgram, "pointLightIntensityMultiplier",
		pointLightIntensityMultiplier);
	labhelper::setUniformSlow(currentShaderProgram, "lightPosition", lightPosition);

	// Simulation parameters
	labhelper::setUniformSlow(currentShaderProgram, "uTime", currentTime);
	labhelper::setUniformSlow(currentShaderProgram, "cloudTime", cloudTime);
	labhelper::setUniformSlow(currentShaderProgram, "planetRadius", planetRadius);
	labhelper::setUniformSlow(currentShaderProgram, "cloudlessDepth", cloudlessDepth);
	labhelper::setUniformSlow(currentShaderProgram, "cloudDepth", cloudDepth);
	labhelper::setUniformSlow(currentShaderProgram, "samplingFalloffDistance", samplingFalloffDistance);
	labhelper::setUniformSlow(currentShaderProgram, "cloudScale", cloudScale);
	labhelper::setUniformSlow(currentShaderProgram, "cloudStepMin", cloudStepMin);
	labhelper::setUniformSlow(currentShaderProgram, "cloudShadowCutoff", cloudShadowCutoff);
	labhelper::setUniformSlow(currentShaderProgram, "cloudShadowIntensity", cloudShadowIntensity);
	labhelper::setUniformSlow(currentShaderProgram, "cloudStepMax", cloudStepMax);
	labhelper::setUniformSlow(currentShaderProgram, "cloudLightingFalloff", cloudLightingFalloff);
	labhelper::setUniformSlow(currentShaderProgram, "atmosphereDepth", atmosphereDepth);
	labhelper::setUniformSlow(currentShaderProgram, "atmosphereDensityFalloff", atmosphereDensityFalloff);
	labhelper::setUniformSlow(currentShaderProgram, "atmosphereDensityAtSeaLevel", atmosphereDensityAtSeaLevel);
	vec3 scatteringCoefficients = vec3(pow(300 / colorBandWavelengths.x, 4), pow(300 / colorBandWavelengths.y, 4), pow(300 / colorBandWavelengths.z, 4)) * atmosphereScatteringStrength;
	labhelper::setUniformSlow(currentShaderProgram, "atmosphereScatteringCoefficients", scatteringCoefficients);
	

	// Sampling
	labhelper::setUniformSlow(currentShaderProgram, "samplingIncreaseFactor", samplingIncreaseFactor);
	labhelper::setUniformSlow(currentShaderProgram, "samplingIncreaseDepth", samplingIncreaseDepth);
	labhelper::setUniformSlow(currentShaderProgram, "samplingFalloffDistance", samplingFalloffDistance);

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

void drawSolidGeometry(GLuint currentShaderProgram, 
	const mat4& viewMatrix,
	const mat4& projectionMatrix,
	const mat4& lightViewMatrix, 
	const mat4& lightProjectionMatrix) {
	glUseProgram(currentShaderProgram);
	glFrontFace(GL_CW); // The models are rendered inside out so we flip what is considered to be the front face

	// Light source
	vec4 viewSpaceLightPosition = viewMatrix * vec4(lightPosition, 1.0f);
	labhelper::setUniformSlow(currentShaderProgram, "point_light_color", pointLightColor);
	labhelper::setUniformSlow(currentShaderProgram, "point_light_intensity_multiplier",
		point_light_intensity_multiplier);
	labhelper::setUniformSlow(currentShaderProgram, "viewSpaceLightPosition", vec3(viewSpaceLightPosition));
	labhelper::setUniformSlow(currentShaderProgram, "viewSpaceLightDir",
		normalize(vec3(viewMatrix * vec4(-lightPosition, 0.0f))));
	labhelper::setUniformSlow(currentShaderProgram, "spotOuterAngle", std::cos(radians(outerSpotlightAngle)));
	labhelper::setUniformSlow(currentShaderProgram, "spotInnerAngle", std::cos(radians(innerSpotlightAngle)));
	labhelper::setUniformSlow(currentShaderProgram, "useSoftFalloff", useSoftFalloff ? 1 : 0);

	glActiveTexture(GL_TEXTURE10);
	glBindTexture(GL_TEXTURE_2D, shadowMapFB.depthBuffer);
	mat4 lightMatrix = translate(vec3(0.5f)) * scale(vec3(0.5f)) * lightProjectionMatrix * lightViewMatrix * inverse(viewMatrix);
	labhelper::setUniformSlow(currentShaderProgram, "lightMatrix", lightMatrix);

	// Planet
	labhelper::setUniformSlow(currentShaderProgram, "modelViewProjectionMatrix",
		projectionMatrix * viewMatrix * planetModelMatrix);
	labhelper::setUniformSlow(currentShaderProgram, "modelViewMatrix", viewMatrix * planetModelMatrix);
	labhelper::setUniformSlow(currentShaderProgram, "normalMatrix",
		inverse(transpose(viewMatrix * planetModelMatrix)));

	labhelper::render(planetModel);
}

void debugDrawLight(const glm::mat4& viewMatrix,
	const glm::mat4& projectionMatrix,
	const glm::vec3& worldSpaceLightPos)
{
	mat4 modelMatrix = glm::translate(worldSpaceLightPos);
	glUseProgram(shaderProgram);
	labhelper::setUniformSlow(shaderProgram, "modelViewProjectionMatrix",
		projectionMatrix * viewMatrix * modelMatrix);
	labhelper::render(planetModel);
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
	//viewMatrix = lookAt(cameraPosition, cameraPosition - cameraDirection, worldUp);

	vec4 lightStartPosition = vec4(40.0f, 40.0f, 0.0f, 1.0f);
	if(animateLight)
		lightPosition = vec3(rotate(currentTime, worldUp) * lightStartPosition);
	mat4 lightViewMatrix = lookAt(lightPosition, vec3(0.0f), worldUp);
	mat4 lightProjMatrix = perspective(radians(45.0f), 1.0f, 25.0f, 100.0f);

	///////////////////////////////////////////////////////////////////////////
	// Set Up Shadow Map
	///////////////////////////////////////////////////////////////////////////
	if (shadowMapFB.width != shadowMapResolution || shadowMapFB.height != shadowMapResolution) {
		shadowMapFB.resize(shadowMapResolution, shadowMapResolution);
	}

	if (shadowMapClampMode == ClampMode::Edge) {
		glBindTexture(GL_TEXTURE_2D, shadowMapFB.depthBuffer);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	}

	if (shadowMapClampMode == ClampMode::Border) {
		glBindTexture(GL_TEXTURE_2D, shadowMapFB.depthBuffer);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
		vec4 border(shadowMapClampBorderShadowed ? 0.f : 1.f);
		glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, &border.x);
	}

	if (useHardwarePCF) {
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	}
	else {
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
	if (usePolygonOffset) {
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(polygonOffset_factor, polygonOffset_units);
	}

	glBindTexture(GL_TEXTURE_2D, shadowMapFB.depthBuffer);
	glBindFramebuffer(GL_FRAMEBUFFER, shadowMapFB.framebufferId);
	glViewport(0, 0, shadowMapFB.width, shadowMapFB.height);
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClearDepth(1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	drawSolidGeometry(depthProgram, viewMatrix, projMatrix, lightViewMatrix, lightProjMatrix);
	if (usePolygonOffset) {
		glDisable(GL_POLYGON_OFFSET_FILL);
	}

	// Draw to fbo from camera
	glBindFramebuffer(GL_FRAMEBUFFER, rasterizedFBO.framebufferId);
	glViewport(0, 0, rasterizedFBO.width, rasterizedFBO.height);
	glClearColor(0, 0, 0, 0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	drawSolidGeometry(shaderProgram, viewMatrix, projMatrix, lightViewMatrix, lightProjMatrix);
	debugDrawLight(viewMatrix, projMatrix, vec3(lightPosition));

	{
		labhelper::perf::Scope s("Ray Marching");
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

		drawScene(raymarchingProgram, viewMatrix, projMatrix);
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
		if (event.type == SDL_MOUSEBUTTONUP && event.button.button == SDL_BUTTON_LEFT && (!labhelper::isGUIvisible() || !ImGui::GetIO().WantCaptureMouse)) {
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
	const uint8_t* state = SDL_GetKeyboardState(nullptr);

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
	ImGui::Text("Simulation");
	ImGui::SliderFloat("Cloud Speed", &cloudMovementSpeed, 0.0f, 1.0f);
	ImGui::SliderFloat("Sampling Factor", &samplingIncreaseFactor, 1.0f, 100.0f);
	ImGui::SliderFloat("Sampling Depth", &samplingIncreaseDepth, 0.0f, 100.0f);
	ImGui::SliderFloat("Sampling Falloff", &samplingFalloffDistance, 0.1f, 100.0f);
	ImGui::SliderFloat("Planet Radius", &planetRadius, 1, 100.0f);
	ImGui::SliderFloat("Cloudless Depth", &cloudlessDepth, 0.1f, 10);
	ImGui::SliderFloat("Cloud Depth", &cloudDepth, 0.0f, 1);
	ImGui::SliderFloat("Cloud Scale", &cloudScale, 0.1f, 1);
	ImGui::SliderFloat("Cloud Step Min", &cloudStepMin, 0, 1);
	ImGui::SliderFloat("Cloud Step Max", &cloudStepMax, 0, 1);
	ImGui::SliderFloat("Cloud Shadow Cutoff", &cloudShadowCutoff, 0, 1);
	ImGui::SliderFloat("Cloud Shadow Intensity", &cloudShadowIntensity, 0, 10);
	ImGui::SliderFloat("Cloud Lighting Fraction", &cloudLightingFalloff, 0.001, 1);
	ImGui::Text("Atmosphere");
	ImGui::SliderFloat("Depth", &atmosphereDepth, 0, 10);
	ImGui::SliderFloat("Density Falloff", &atmosphereDensityFalloff, 0, 10);
	ImGui::SliderFloat("Density at Sea Level", &atmosphereDensityAtSeaLevel, 0, 1);
	ImGui::SliderFloat3("Scattering Wavelengths", (float*)&colorBandWavelengths, 0, 1000);
	ImGui::SliderFloat("Scattering Strength", &atmosphereScatteringStrength, 0, 10);
	ImGui::Checkbox("Animate light", &animateLight);
	ImGui::SliderFloat3("Sun Position", (float*)&lightPosition, -200, 200);
	ImGui::Text("Controls");
	ImGui::SliderFloat("Camera Speed", &cameraSpeed, 0.1f, 100.0f);
	labhelper::perf::drawEventsWindow();
}

int main(int argc, char* argv[])
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

	labhelper::freeModel(planetModel);

	// Shut down everything. This includes the window and all other subsystems.
	labhelper::shutDown(g_window);
	return 0;
}