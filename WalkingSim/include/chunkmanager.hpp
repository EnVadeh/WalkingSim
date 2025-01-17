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
	GLint length = 10;
	GLint breadth = 10;
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


vector2D randomGradient(vector2D p) {
	p = p + 0.1;
	float x = dot(p, vector2D(123.4, 234.5));
	float y = dot(p, vector2D(234.5, 345.6));
	vector2D gradient = vector2D(x, y);
	gradient.sinof();
	gradient = gradient * 43758.5453;

	gradient.sinof();
	return gradient;
}

template<typename T>
void perlinNoise(image2D<T>& image) {
	float size_x = image.size.x;
	float size_y = image.size.y;
	image2D<vector2D> PN(size_x, size_y); //For each pixel, we need 4 constantish values right? So x,y is itself, the corner to it's right is x+1,y value, the corner to its up/down depending on the coordinate system is it's y+1/y-1, and the corner to its diagonal is x+1,y+1
	for (size_t x; x < size_x; x++) {
		for (size_t y; y < size_t; y++) {
			PN.write(x, y, randomGradient(vector2D(x, y)));
		}
	}
}