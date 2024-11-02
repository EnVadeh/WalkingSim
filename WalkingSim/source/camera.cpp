#include "camera.hpp"

camera::camera(glm::vec3 initPos, glm::vec3 initDir, GLfloat yaw, GLfloat pitch) : vEye(initPos), vFront(initDir), fYaw(yaw), fPitch(pitch) {
	matView = glm::lookAt(vEye, vFront, vUp);
	matProj = glm::perspective(float(glm::radians(60.0f)), 1.0f, 0.1f, 10000.0f); //could make my own frustum but this is far simpler and more effective afaik
	setCam();
}

void camera::setCam() {
	matProjView = matProj * matView;

	//std::cout << "The position is: (" << vEye.x << "," << vEye.y << ","<< vEye.z << ")" << std::endl;
	//std::cout << "The front is: (" << vFront.x << "," << vFront.y << ","<< vFront.z << ")" << std::endl;
	//std::cout << "The up is: (" << vUp.x << "," << vUp.y << ","<< vUp.z << ")" << std::endl;
	for (size_t i = 0; i < SHADER_COUNT; i++) {
		setUniform(SHADERS[i], "matProjView", matProjView);
		setUniform(SHADERS[i], "matProj", matProj);
		setUniform(SHADERS[i], "matView", matView);
		setUniform(SHADERS[i], "vCamPos", vEye);
	}
}

void camera::processKeyboardInput(GLint key) {
	float cameraSpeed = 7.5f * deltaTime;
	if (key == GLFW_KEY_W) vEye += cameraSpeed * vFront;
	if (key == GLFW_KEY_S) vEye -= cameraSpeed * vFront;
	if (key == GLFW_KEY_A) vEye -= glm::normalize(glm::cross(vFront, vUp)) * cameraSpeed;
	if (key == GLFW_KEY_D) vEye += glm::normalize(glm::cross(vFront, vUp)) * cameraSpeed;
	matView = glm::lookAt(vEye, vFront + vEye, vUp);
	setCam();
}

void camera::processMouseInput(GLfloat xoffset, GLfloat yoffset) {
	float sensitivity = 0.1f;
	xoffset *= sensitivity;
	yoffset *= sensitivity;

	fYaw += xoffset;
	fPitch += yoffset;

	if (fPitch > 89.0f) 
		fPitch = 89.0f;
	if (fPitch < -89.0f) 
		fPitch = -89.0f;

	updateCameraVectors();
	matView = glm::lookAt(vEye, vFront + vEye, vUp);
	setCam();
}

void camera::updateCameraVectors() {
	glm::vec3 frontVec;
	frontVec.x = cos(glm::radians(fYaw)) * cos(glm::radians(fPitch));
	frontVec.y = sin(glm::radians(fPitch));
	frontVec.z = sin(glm::radians(fYaw)) * cos(glm::radians(fPitch));
	vFront = glm::normalize(frontVec);
	glm::vec3 vRight = glm::normalize(glm::cross(vFront, glm::vec3(0.0f, 1.0f, 0.0f)));
	vUp = glm::normalize(glm::cross(vRight, vFront));
}
