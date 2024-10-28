#pragma once

#include "utils.hpp"
#include "buffer.hpp"
#include "camera.hpp"

class chunkManager { //I think I should save some chunks and draw them, and discard them later if I go too far away
private:
	GLint vPos[2];
	std::vector<glm::vec3> chunks;
	std::shared_ptr<terrain> chunk;
	camera* cCam = nullptr;
	bool negative_x = false;
	bool negative_z = false;
	GLint length = 10;
	GLint breadth = 10;
public:
	chunkManager(camera &cam);
	void checkPos(GLuint shaderID);
	void draw(GLuint shaderID);
};