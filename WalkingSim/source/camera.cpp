#include "camera.hpp"

camera::camera(glm::vec3 initPos, glm::vec3 initDir) {
	vCenter = initDir;
	vEye = initPos;
	matView = glm::lookAt(vEye, vCenter, vUp);
	matProj = glm::perspective(float(glm::radians(60.0f)), 1.0f, 0.1f, 10000.0f); //could make my own frustum but this is far simpler and more effective afaik
	matProjView = matProj * matView;	
	setCam();
}

void camera::setCam() {
	for (size_t i = 0; i < SHADER_COUNT; i++)
		setUniform(SHADERS[i], "matProjView", matProjView);
}