#include "chunkmanager.hpp"

chunkManager::chunkManager(camera& cam) {
	cCam = &cam;
	chunk = std::make_shared<terrain>(length, breadth);
}

void chunkManager::checkPos(GLuint shaderID) {

	vPos[0] = cCam->vEye.x / (length);
	vPos[1] = cCam->vEye.z / (breadth);

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

	//Center Chunk's coordinates.
	tempV.x = chunks[0].x; 
	tempV.z = chunks[0].z;

	//if moving towards negative x or negative z axis then the chunk's coordinates is decreased
	if (negative_x == true)
		tempV.x = tempV.x - 1;

	if (negative_z == true) 
		tempV.z = tempV.z - 1;
	
	//Scale it up to the chunk size
	tempV.x = tempV.x * length;
	tempV.y = 0;
	tempV.z = tempV.z * breadth;

	glm::vec3 original; //named original because this is the new point that other temporary positions will be relative to
	original = tempV;

	//Doing a loop where 9 chunks are rendered 
	for (int i = 0; i < 3; i++) {
		tempV.z = (i - 1) * breadth + original.z; //row wise, offset from the new z's position
		for (int j = 0; j < 3; j++) {
			tempV.x = (j - 1) * length + original.x; //column wise, offset form the new x's position
			chunk->draw(shaderID, tempV, glm::vec3(1, 1, 1)); //draw each time offset from the center (original) position
		}
	}
}

atmosphereLUTs::atmosphereLUTs(const atmosphereParams& params) : atmosphere(params) {
	initializeLUTs();
}

void atmosphereLUTs::initializeLUTs() {
	std::vector<atmosphereParams> temp;
	temp.push_back(atmosphere);
	uniformBuffer<atmosphereParams> atmosphereBuff(temp, drawFreq::staticDraw);
	atmosphereBuff.bind();
	const float EarthRayleighScaleHeight = 8.0f;
	const float EarthMieScaleHeight = 1.2f;
	DP.push_back({ 25000.0, 0.0, 0.0, 1.0 / 15000.0, -2.0 / 3.0});
	DP.push_back({ 0.0, 0.0, 0.0, -1.0 / 15000.0, 8.0 / 3.0 });
	DP.push_back({ 0.0, 1.0, -1.0/ 8000.0f, 0, 0.0 });//rayleigh //in rayleigh densityprifle, 0 is the 0 then 1 is this
	DP.push_back({ 0.0, 1.0, -1.0/ 1200.0f, 0, 0.0 });//mie
	DP.push_back({ 0.0, 0.0, 0.0 , 0, 0.0 });//empty
	

	uniformBuffer<densityProfileLayer> profileBuff(DP, drawFreq::staticDraw);
	profileBuff.bind();
	createTransmittanceLUT();
	createScatteringLUT();
	lutNames.push_back(a1);
	lutNames.push_back(a2);
	lutNames.push_back(a3);
	lutNames.push_back(a4);
	lutNames.push_back(a5);
	lutNames.push_back(a6);
	lutNames.push_back(a7);
	lutNames.push_back(a8);
}

void atmosphereLUTs::bind(GLuint shaderID) {
	for (u8 i = 0; i < 8; i++)
		CO[i].bind(shaderID, i, lutNames[i]);
}

void atmosphereLUTs::bind(GLuint shaderID, GLuint from, GLuint to) {
	for (u8 i = from; i < to; i++)
		CO[i].bind(shaderID, i, lutNames[i]);
}

void atmosphereLUTs::createTransmittanceLUT() {
	//CO[0].setup(TRANSMITTANCE_W, TRANSMITTANCE_H, 0);
	CO[0].setup(1024, 1024, 0);
}

void atmosphereLUTs::createScatteringLUT() {

	for (u8 i = 1; i < 4; i++)
		CO[i].setup(SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT, SCATTERING_TEXTURE_DEPTH, i); // 1 = scatteringLUT, 2 = rayleighDeltaLUT, 3 = mieDeltaLUT
	CO[4].setup(IRRADIANCE_TEXTURE_WIDTH, IRRADIANCE_TEXTURE_HEIGHT, 4); // 4 = directIrradiance
	CO[5].setup(SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT, SCATTERING_TEXTURE_DEPTH, 5); // 5 = scatteringDesntiyLUT
	CO[6].setup(IRRADIANCE_TEXTURE_WIDTH, IRRADIANCE_TEXTURE_HEIGHT, 6); // 6 = indirectIrradiance
	CO[7].setup(SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT, SCATTERING_TEXTURE_DEPTH, 7); // 7 = multiScatteringLUT

}
