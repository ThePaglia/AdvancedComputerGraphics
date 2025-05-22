#version 420

// required by GLSL spec Sect 4.5.3 (though nvidia does not, amd does)
precision highp float;

///////////////////////////////////////////////////////////////////////////////
// Material
///////////////////////////////////////////////////////////////////////////////
uniform vec3 material_color;
uniform vec3 material_reflectivity;
uniform float material_metalness;
uniform float material_fresnel;
uniform float material_shininess;
uniform vec3 material_emission;

uniform int has_color_texture;
layout(binding = 0) uniform sampler2D colorMap;
uniform int has_emission_texture;
layout(binding = 5) uniform sampler2D emissiveMap;

///////////////////////////////////////////////////////////////////////////////
// Directional light source (Sun)
///////////////////////////////////////////////////////////////////////////////
uniform vec3 directional_light_color = vec3(1.0, 1.0, 1.0);
uniform float directional_light_intensity_multiplier = 50.0;

///////////////////////////////////////////////////////////////////////////////
// Constants
///////////////////////////////////////////////////////////////////////////////
#define PI 3.14159265359

///////////////////////////////////////////////////////////////////////////////
// Input varyings from vertex shader
///////////////////////////////////////////////////////////////////////////////
in vec2 texCoord;
in vec3 viewSpaceNormal;
in vec3 viewSpacePosition;

///////////////////////////////////////////////////////////////////////////////
// Input uniform variables
///////////////////////////////////////////////////////////////////////////////
uniform mat4 viewInverse;
uniform vec3 viewSpaceLightPosition;
uniform vec3 planetOrigin = vec3(0, 0, 0);

///////////////////////////////////////////////////////////////////////////////
// Output color
///////////////////////////////////////////////////////////////////////////////
layout(location = 0) out vec4 fragmentColor;

// Shadows
in vec4 shadowMapCoord;
layout(binding = 10) uniform sampler2DShadow shadowMapTex;
uniform vec3 viewSpaceLightDir;


vec3 calculateDirectIllumiunation(vec3 wo, vec3 n, vec3 base_color)
{
	vec3 direct_illum = base_color;

	// Directional light term (the sun), distance is not relevant here
	vec3 Li = directional_light_intensity_multiplier * directional_light_color;
	vec3 wi = normalize(viewSpaceLightPosition - viewSpacePosition);
	if(dot(wi, n) <= 0.0)
		return vec3(0.0);
	
	vec3 diffuse_term = base_color * (1.0 / PI) * dot(n, wi) * Li;
	direct_illum = diffuse_term;

	// Apply gamma correction
	return pow(direct_illum, vec3(1.0 / 2.2));
}

vec3 calculateIndirectIllumination(vec3 wo, vec3 n, vec3 base_color)
{
	// The idea here is to create a vector from the planet's origin towards the fragment and compare that vector to the sun direction
	// This can be done in world space, but to save precious computation, we do it in view space instead (as we have those variables already)
	// The planet's origin in view space. Since we get the other ones in view space, this is the only one we need to transform!
	vec3 viewSpacePlanetOrigin = vec3(inverse(viewInverse) * vec4(planetOrigin, 1));
	// The direction from the planet's center towards this fragment in view space
	vec3 planetNormal = normalize(viewSpacePosition - viewSpacePlanetOrigin);

	// The indirect illumination (i.e. shadow strength) is scaled by where the fragment is on the planet (or actually where it would be on a perfect sphere)
	// This acts similarly to lighting a perfect sphere diffusely. The negative sign is needed since we compare towards the sun
	float diffuse = max(dot(planetNormal, -viewSpaceLightDir), 0);

	return base_color * diffuse;
}

void main()
{
	float visibility = 1.0;
	float attenuation = 1.0;

	// TODO: Figure out why visibility is not correctly calculated
	visibility = textureProj(shadowMapTex, shadowMapCoord);

	vec3 wo = -normalize(viewSpacePosition);
	vec3 n = normalize(viewSpaceNormal);

	vec3 base_color = material_color;
	if(has_color_texture == 1)
	{
		base_color = texture(colorMap, texCoord).rgb;
	}

	// Direct illumination
	vec3 direct_illumination_term = visibility * calculateDirectIllumiunation(wo, n, base_color);

	// Indirect illumination
	vec3 indirect_illumination_term = calculateIndirectIllumination(wo, n, base_color);

	///////////////////////////////////////////////////////////////////////////
	// Add emissive term. If emissive texture exists, sample this term.
	///////////////////////////////////////////////////////////////////////////
	vec3 emission_term = material_emission;
	if(has_emission_texture == 1)
	{
		emission_term = texture(emissiveMap, texCoord).rgb;
	}

	vec3 shading = direct_illumination_term + indirect_illumination_term + emission_term;

	fragmentColor = vec4(shading, 1.0);
	return;
}
