#pragma once

#include "utils.hpp"
#include "buffer.hpp"
#include "camera.hpp"

class chunkManager {
private:
	GLint vPos[2];
	std::vector<glm::vec2> chunks;
	std::shared_ptr<terrain> chunk;
	camera* cCam = nullptr;
public:
	chunkManager(camera &cam);
	void checkPos();
	void checkk();
};