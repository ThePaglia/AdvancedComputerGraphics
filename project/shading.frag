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
// Environment
///////////////////////////////////////////////////////////////////////////////
layout(binding = 6) uniform sampler2D environmentMap;
layout(binding = 7) uniform sampler2D irradianceMap;
layout(binding = 8) uniform sampler2D reflectionMap;
uniform float environment_multiplier;

///////////////////////////////////////////////////////////////////////////////
// Light source
///////////////////////////////////////////////////////////////////////////////
uniform vec3 point_light_color = vec3(1.0, 1.0, 1.0);
uniform float point_light_intensity_multiplier = 50.0;

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

///////////////////////////////////////////////////////////////////////////////
// Output color
///////////////////////////////////////////////////////////////////////////////
layout(location = 0) out vec4 fragmentColor;

// Shadows
in vec4 shadowMapCoord;
layout(binding = 10) uniform sampler2DShadow shadowMapTex;
uniform vec3 viewSpaceLightDir;
uniform float spotOuterAngle;
uniform float spotInnerAngle;
uniform int useSoftFalloff;


vec3 calculateDirectIllumiunation(vec3 wo, vec3 n, vec3 base_color)
{
	vec3 direct_illum = base_color;
	///////////////////////////////////////////////////////////////////////////
	// Task 1.2 - Calculate the radiance Li from the light, and the direction
	//            to the light. If the light is backfacing the triangle,
	//            return vec3(0);
	///////////////////////////////////////////////////////////////////////////
	const float distance_to_light = length(viewSpaceLightPosition - viewSpacePosition);
	const float falloff_factor = 1.0 / (distance_to_light * distance_to_light);
	vec3 Li = point_light_intensity_multiplier * point_light_color * falloff_factor;
	vec3 wi = normalize(viewSpaceLightPosition - viewSpacePosition);
	if(dot(wi, n) <= 0.0)
		return vec3(0.0);

	///////////////////////////////////////////////////////////////////////////
	// Task 1.3 - Calculate the diffuse term and return that as the result
	///////////////////////////////////////////////////////////////////////////
	
	vec3 diffuse_term = base_color * (1.0 / PI) * dot(n, wi) * Li;
	direct_illum = diffuse_term;

	return pow(direct_illum, vec3(1.0 / 2.2));
}

vec3 calculateIndirectIllumination(vec3 wo, vec3 n, vec3 base_color)
{
	// TODO: Add indirect lighting according to the planet's normal at this point
	return vec3(0.0);
}

void main()
{
	float visibility = 1.0;
	float attenuation = 1.0;

	// TODO: Figure out why visibility is not correctly calculated
	visibility = textureProj(shadowMapTex, shadowMapCoord);
		
	if(useSoftFalloff == 1) {
		vec3 posToLight = normalize(viewSpaceLightPosition - viewSpacePosition);
		float cosAngle = dot(posToLight, -viewSpaceLightDir);

		// Spotlight with hard border:
		float spotAttenuation = smoothstep(spotOuterAngle, spotInnerAngle, cosAngle);
		visibility *= spotAttenuation;
	}

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
