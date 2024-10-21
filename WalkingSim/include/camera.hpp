#pragma once
#include "utils.hpp"

class camera {
private:
	glm::vec3 vEye;
	glm::vec3 vUp = {0.0f, 1.0f, 0.0f};
	glm::vec3 vCenter;
	glm::mat4 matView;
	glm::mat4 matProj;
	glm::mat4 matProjView; 
public:
	camera(glm::vec3 initPos, glm::vec3 initDir);
	void keyMovement();
	void mouseMovement();
	void setCam();
};