#pragma once
#include "utils.hpp"

class camera {
private:
	glm::vec3 vEye;
	glm::vec3 vUp = {0.0f, 1.0f, 0.0f};
	glm::vec3 vFront;
	glm::mat4 matView;
	glm::mat4 matProj;
	glm::mat4 matProjView;

	void updateCameraVectors();
public:
	GLfloat fYaw = 0; //horizontal
	GLfloat fPitch = 0; //vertical
	
	camera(glm::vec3 initPos, glm::vec3 initDir, GLfloat yaw, GLfloat pitch);
	void processKeyboardInput(GLint key);
	void processMouseInput(GLfloat xoffset, GLfloat yoffset);
	void setCam();
};