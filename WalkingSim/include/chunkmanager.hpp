#pragma once
#include "utils.hpp"
#include "buffer.hpp"
#include "camera.hpp"
#include "texture.hpp"


class chunkManager { //I think I should save some chunks and draw them, and discard them later if I go too far away
private:
	GLint vPos[2];
	std::vector<glm::vec3> chunks;
	std::shared_ptr<terrain> chunk;
	camera* cCam = nullptr;
	bool negative_x = false;
	bool negative_z = false;
	GLint length = 20;
	GLint breadth = 20;
	void draw(GLuint shaderID);
public:
	chunkManager(camera &cam);
	void checkPos(GLuint shaderID);
};

struct atmosphereParams { //Precomputed Atmospheric Scattering Eric Bruneton, Fabrice Neyret
	float earthRad = 6360e3f;//6371e3f;
	float atmosphereRad = 6420e3f;//6471e3f; //100 km atmosphere
	float Hr = 8000.0f; //Rayleight scale height //How the density of air molecules with height
	float Hm = 1200.0f; //Mei scale height //How the density of aerosols scales with height
	//--16 byte alignment--//
	//For RGB wavelengths = 612, 549, 465 nm (self calculated)
	//glm::vec3 betaR = {8.84e-6f, 1.365e-5f, 2.65e-5f};
	glm::vec3 betaR = {5.84e-6f, 1.365e-5f, 3.31e-5f};
	float padding1;
	//--16 byte alignment--//
	//glm::vec3 betaM = {2.1e-5f, 2.1e-5f, 2.1e-5f};
	glm::vec3 betaM = {4.0e-5f, 4.0e-5f, 4.0e-5f};
	float padding2;
	//--16 bit alignment--//
	float meiG = 0.76f;
	float padding3;
	float padding4;
	float padding5;
	//--16 bit alignment--//
	//glm::vec3 betaMext = { 2.33e-5, 2.33e-5, 2.33e-5 };
	glm::vec3 betaMext = { 2.33e-5, 2.33e-5, 2.33e-5 };
	float padding6;
	//--16 bit alignment--//
};

class atmosphereLUTs {
private:
	std::vector<densityProfileLayer> DP;
	atmosphereParams atmosphere;
	computeOutput CO[8];
	std::vector<std::string>lutNames;
	std::string a1 = "transmittanceLUT"; 
	std::string a2 = "scatteringLUT"; 
	std::string a3 = "rayleighLUT";
	std::string a4 = "mieLUT";
	std::string a5 = "deltaIrradianceLUT"; //as eric bruneton paper: 16 x 64 //cehck again not compeltely sure
	std::string a6 = "scatteringDensityLUT";
	std::string a7 = "irradianceLUT";
	std::string a8 = "multiScatteringLUT";

	static const int TRANSMITTANCE_W = 256; 
	static const int TRANSMITTANCE_H = 64;
	static const int SCATTERING_R = 32;
	static const int SCATTERING_MU = 128;
	static const int SCATTERING_MU_S = 32;
	static const int SCATTERING_NU = 8;
	static const int SCATTERING_TEXTURE_WIDTH = SCATTERING_MU_S * SCATTERING_NU;
	static const int SCATTERING_TEXTURE_HEIGHT = SCATTERING_MU;
	static const int SCATTERING_TEXTURE_DEPTH = SCATTERING_R;
	static const int IRRADIANCE_TEXTURE_WIDTH = 64;
	static const int IRRADIANCE_TEXTURE_HEIGHT = 16;

	void initializeLUTs();
	void createTransmittanceLUT();
	void createScatteringLUT();

public:
	atmosphereLUTs(const atmosphereParams& params);
	void bind(GLuint shaderID);
	void bind(GLuint shaderID, GLuint from, GLuint to);
};

