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
ivec2 g_prevMouseCoords = { -1, -1 };
bool g_isMouseDragging = false;

// Shader programs
GLuint raymarchingProgram;

// Camera parameters.
vec3 worldUp(0.0f, 1.0f, 0.0f);
vec3 cameraPosition(0.0f, 0.0f, 5.0f);
vec3 cameraDirection = normalize(vec3(0.0f) - cameraPosition);
vec3 cameraRight = cross(cameraDirection, worldUp);
vec3 cameraUp = cross(cameraRight, cameraDirection);

float cameraSpeed = 10.f;

// Texture parameters
GLuint noiseTexture;

// Light parameters
vec3 lightPosition;
vec3 pointLightColor = vec3(1.f, 1.f, 1.f);

float pointLightIntensityMultiplier = 10000.0f;

void loadShaders(bool is_reload)
{
	GLuint shader = labhelper::loadShaderProgram("../project/raymarching.vert", "../project/raymarching.frag", is_reload);
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

// This function is called once at the start of the program and never again
void initialize()
{
	ENSURE_INITIALIZE_ONLY_ONCE();

	// Load Shaders
	loadShaders(false);

	// Load noise texture
	loadNoiseTexture("../textures/noise.png");

	glEnable(GL_DEPTH_TEST); // enable Z-buffering
	glEnable(GL_CULL_FACE);	 // enables backface culling
}

// This function is used to draw the main objects on the scene
void drawScene(GLuint currentShaderProgram,
	const mat4& viewMatrix,
	const mat4& projectionMatrix)
{
	glUseProgram(currentShaderProgram);

	// Bind the noise texture
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, noiseTexture);
	labhelper::setUniformSlow(currentShaderProgram, "uNoiseTexture", 0);

	// Light source
	vec4 viewSpaceLightPosition = viewMatrix * vec4(lightPosition, 1.0f);
	labhelper::setUniformSlow(currentShaderProgram, "pointLightColor", pointLightColor);
	labhelper::setUniformSlow(currentShaderProgram, "pointLightIntensityMultiplier",
		pointLightIntensityMultiplier);
	labhelper::setUniformSlow(currentShaderProgram, "viewSpaceLightPosition", vec3(viewSpaceLightPosition));
	labhelper::setUniformSlow(currentShaderProgram, "viewSpaceLightDir",
		normalize(vec3(viewMatrix * vec4(-lightPosition, 0.0f))));

	// uTime
	labhelper::setUniformSlow(currentShaderProgram, "uTime", currentTime);

	// uResolution
	labhelper::setUniformSlow(currentShaderProgram, "uResolution", vec2(windowWidth, windowHeight));

	// Camera
	labhelper::setUniformSlow(currentShaderProgram, "uCameraPos", cameraPosition);
	labhelper::setUniformSlow(currentShaderProgram, "uCameraDir", cameraDirection);
	labhelper::setUniformSlow(currentShaderProgram, "uCameraUp", cameraUp);
	labhelper::setUniformSlow(currentShaderProgram, "uCameraRight", cameraRight);

	labhelper::drawFullScreenQuad();
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
		}
	}

	// setup matrices
	mat4 projMatrix = perspective(radians(45.0f), float(windowWidth) / float(windowHeight), 5.0f, 2000.0f);
	mat4 viewMatrix = lookAt(cameraPosition, cameraPosition + cameraDirection, worldUp);

	vec4 lightStartPosition = vec4(0.0f, 0.0f, 0.0f, 1.0f);
	lightPosition = vec3(rotate(currentTime, worldUp) * lightStartPosition);

	// Draw from camera
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glViewport(0, 0, windowWidth, windowHeight);
	glClearColor(0.2f, 0.2f, 0.8f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Draw the scene
	drawScene(raymarchingProgram, viewMatrix, projMatrix);
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
			g_isMouseDragging = true;
			int x;
			int y;
			SDL_GetMouseState(&x, &y);
			g_prevMouseCoords.x = x;
			g_prevMouseCoords.y = y;
		}

		if (!(SDL_GetMouseState(NULL, NULL) & SDL_BUTTON(SDL_BUTTON_LEFT)))
		{
			g_isMouseDragging = false;
		}

		if (event.type == SDL_MOUSEMOTION && g_isMouseDragging)
		{
			// More info at https://wiki.libsdl.org/SDL_MouseMotionEvent
			int delta_x = event.motion.x - g_prevMouseCoords.x;
			int delta_y = event.motion.y - g_prevMouseCoords.y;
			float rotationSpeed = 0.1f;
			mat4 yaw = rotate(rotationSpeed * deltaTime * -delta_x, worldUp);
			cameraRight = normalize(cross(cameraDirection, worldUp));
			mat4 pitch = rotate(rotationSpeed * deltaTime * -delta_y,
				cameraRight);
			cameraDirection = vec3(yaw * vec4(cameraDirection, 0.0f));
			g_prevMouseCoords.x = event.motion.x;
			g_prevMouseCoords.y = event.motion.y;
		}
	}

	// check keyboard state (which keys are still pressed)
	const uint8_t* state = SDL_GetKeyboardState(nullptr);

	if (state[SDL_SCANCODE_W])
	{
		cameraPosition += cameraSpeed * deltaTime * cameraDirection;
	}
	if (state[SDL_SCANCODE_S])
	{
		cameraPosition -= cameraSpeed * deltaTime * cameraDirection;
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
		//cameraPosition += cameraSpeed * deltaTime * normalize(cross(cameraDirection, worldUp));
	}
	if (state[SDL_SCANCODE_Q])
	{
		//cameraPosition -= cameraSpeed * deltaTime * normalize(cross(cameraDirection, worldUp));
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