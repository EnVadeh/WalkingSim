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
	void draw(GLuint shaderID);
public:
	chunkManager(camera &cam);
	void checkPos(GLuint shaderID);
};

struct atmosphereParams { //Precomputed Atmospheric Scattering Eric Bruneton, Fabrice Neyret
	float earthRad = 6371e3f;
	float atmosphereRad = 6471e3f; //100 km atmosphere
	float Hr = 8000.0f; //Rayleight scale height //How the density of air molecules with height
	float Hm = 1200.0f; //Mei scale height //How the density of aerosols scales with height
	
	//For RGB wavelengths = 680, 550, 440 nm
	glm::vec3 betaR = {5.8e-6f, 13.5e-6f, 331.e-6f};
	glm::vec3 betaM = {2.1e-5f, 2.1e-5f, 2.1e-5f};
	float meiG = 0.76f;
};

class atmosphereLUTs {
private:
	atmosphereParams atmosphere;
	GLuint transmittenceLUT; //as eric bruneton paper: 64 x 256
	GLuint irradianceLUT; //as eric bruneton paper: 16 x 64
	GLuint scatteringLUT; //as eric bruneton paper: 32 x 128 x 32 x 8 



	static const int TRANSMITTANCE_W = 32; 
	static const int TRANSMITTANCE_H = 32;
	static const int SCATTERING_R = 16;
	static const int SCATTERING_MU = 32;
	static const int SCATTERING_MU_S = 16;
	static const int SCATTERING_NU = 4;

	void initializeLUTs();
	void createTransmittanceLUT();
	void createScatteringLUT();
	void uvtoTransmittanceParams(float u, float v, float& r, float& mu);

	glm::vec3 computeTransmittance(float r, float mu);
	glm::vec3 computeScattering(float r, float mu, float mu_s, float nu);

	float rayIntersectSphere(float r, float mu, float radius);

public:
	atmosphereLUTs(const atmosphereParams& params);
	void bind(GLuint shaderID);
};