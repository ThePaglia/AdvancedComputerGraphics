#include <glm/glm.hpp>

namespace procedural
{
	float getHeightOnUnitSphere(glm::vec3 point) {
		return 1.0f - sin(point.x * 10) * 0.1f;
	}
}