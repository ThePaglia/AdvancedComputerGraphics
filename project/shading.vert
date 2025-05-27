#version 420
///////////////////////////////////////////////////////////////////////////////
// Input vertex attributes
///////////////////////////////////////////////////////////////////////////////
layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normalIn;
layout(location = 2) in vec2 texCoordIn;
layout(location = 3) in vec3 colorIn;

///////////////////////////////////////////////////////////////////////////////
// Input uniform variables
///////////////////////////////////////////////////////////////////////////////
uniform mat4 normalMatrix;
uniform mat4 modelViewMatrix;
uniform mat4 modelViewProjectionMatrix;
uniform mat4 lightMatrix;

///////////////////////////////////////////////////////////////////////////////
// Output to fragment shader
///////////////////////////////////////////////////////////////////////////////
out vec2 texCoord;
out vec3 viewSpaceNormal;
out vec3 viewSpacePosition;
out vec4 shadowMapCoord;
out vec3 vertexColor;

void main() {
	gl_Position = modelViewProjectionMatrix * vec4(position, 1.0);
	texCoord = texCoordIn;
	viewSpaceNormal = (normalMatrix * vec4(normalIn, 0.0)).xyz;
	viewSpacePosition = (modelViewMatrix * vec4(position, 1.0)).xyz;
	shadowMapCoord = lightMatrix * vec4(viewSpacePosition, 1.0);
	vertexColor = colorIn;
}
