#include "chunkmanager.hpp"

chunkManager::chunkManager(camera& cam) {
	cCam = &cam;
	chunk = std::make_shared<terrain>(10, 10);
}

void chunkManager::checkPos(GLuint shaderID) { //Entire position thing is fucked man
	vPos[0] = cCam->vEye.x / 10;
	vPos[1] = cCam->vEye.z / 10;

	glm::vec3 temp = { vPos[0], 0, vPos[1]};
	chunks.clear();
	chunks.push_back(temp);
	draw(shaderID);
}

void chunkManager::draw(GLuint shaderID) {
	glm::vec3 tempV;
	tempV.x = chunks[0].x * 10;
	tempV.y = 0;
	tempV.z = chunks[0].z * 10;
	chunk->draw(shaderID, tempV, glm::vec3(1, 1, 1));
}