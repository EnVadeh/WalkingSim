#include "chunkmanager.hpp"

chunkManager::chunkManager(camera& cam) {
	cCam = &cam;
	chunk = std::make_shared<terrain>(length, breadth);
}

void chunkManager::checkPos(GLuint shaderID) {

	vPos[0] = cCam->vEye.x / length;
	vPos[1] = cCam->vEye.z / breadth;

	glm::vec3 temp = { vPos[0], 0, vPos[1]};
	 
	negative_x = false;
	negative_z = false;
	
	if (cCam->vEye.x < 0)
		negative_x = true;
	if (cCam->vEye.z < 0)
		negative_z = true;

	chunks.clear();
	chunks.push_back(temp);
	draw(shaderID);
}

void chunkManager::draw(GLuint shaderID) {
	glm::vec3 tempV;
	tempV.x = chunks[0].x;
	tempV.z = chunks[0].z;

	if (negative_x == true)
		tempV.x = tempV.x - 1;

	if (negative_z == true) 
		tempV.z = tempV.z - 1;
	
	tempV.x = tempV.x * length;
	tempV.y = 0;
	tempV.z = tempV.z * breadth;

	glm::vec3 original;
	original = tempV;

	for (int i = 0; i < 3; i++) {
		tempV.z = (i - 1) * breadth + original.z;
		for (int j = 0; j < 3; j++) {
			tempV.x = (j - 1) * length + original.x;
			chunk->draw(shaderID, tempV, glm::vec3(1, 1, 1));
		}
	}

}