#include "chunkmanager.hpp"

chunkManager::chunkManager(camera& cam) {
	cCam = &cam;
	chunks.resize(9);

	//chunk = std::make_shared<terrain>(10, 10);
}

void chunkManager::checkPos() { //Entire position thing is fucked man
	vPos[0] = cCam->vEye.x / 10;
	vPos[1] = cCam->vEye.z / 10;

	glm::vec2 temp;

	temp.x = vPos[0] - 1;
	temp.y = vPos[1] - 1;
	chunks.push_back(temp);
	temp.x = vPos[0];
	temp.y = vPos[1] - 1;
	chunks.push_back(temp);
	temp.x = vPos[0] + 1;
	temp.y = vPos[1] - 1;

	chunks.push_back(temp);
	temp.x = vPos[0] - 1;
	temp.y = vPos[1];
	chunks.push_back(temp);
	temp.x = vPos[0];
	temp.y = vPos[1];
	chunks.push_back(temp);
	temp.x = vPos[0] + 1;
	temp.y = vPos[1];
	chunks.push_back(temp);
	temp.x = vPos[0] - 1;
	temp.y = vPos[1] + 1;
	chunks.push_back(temp);
	temp.x = vPos[0];
	temp.y = vPos[1] + 1;
	chunks.push_back(temp);
	temp.x = vPos[0] + 1;
	temp.y = vPos[1] + 1;
	chunks.push_back(temp);

	//::cout << "Position is from object directly : ( " << cCam->vEye.x << " , " << cCam->vEye.z << " )" << std::endl;
	//std::cout << "Position is: ( " << vPos[0] << " , " << vPos[1] << " )" << std::endl;
	//for (size_t i = 0; i < 9; i++)
	//	std::cout << "The positions are: " << chunks[i].x << " , " << chunks[i].y << std::endl;
	
}

void chunkManager::checkk() {
}