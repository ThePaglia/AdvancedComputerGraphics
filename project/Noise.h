#pragma once
#include <glm/glm.hpp>

namespace procedural
{
	class Noise
	{
	public:
	};

	float getHeightOnUnitSphere(glm::vec3 point);
	glm::vec3 getColorForHeight(float height);
}