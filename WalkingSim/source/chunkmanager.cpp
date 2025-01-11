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
	//Remnant from when I tried precomputing in the CPU FOR SOME REASON
	//std::vector<float> data(TRANSMITTANCE_W * TRANSMITTANCE_H * 3);

	//#pragma omp parallel for collapse(2) //nested fors(2) parallelised
	//for (int j = 0; j < TRANSMITTANCE_H; j++)
	//	for (int i = 0; i < TRANSMITTANCE_W; i++) {
	//		float r, mu;
	//		uvtoTransmittanceParams(float(i) / (TRANSMITTANCE_W - 1), float(j) / (TRANSMITTANCE_H - 1), r, mu);
	//		
	//		glm::vec3 transmittance = computeTransmittance(r, mu);
	//		int idx = (j * TRANSMITTANCE_W + i) * 3;
	//		data[idx + 0] = transmittance.x;
	//		data[idx + 1] = transmittance.y;
	//		data[idx + 2] = transmittance.z;
	//	}
	//CO[0].setup(TRANSMITTANCE_W, TRANSMITTANCE_H, 0);
	CO[0].setup(1024, 1024, 0);
}

void atmosphereLUTs::createScatteringLUT() {
	//std::vector<float> data(SCATTERING_R * SCATTERING_MU * SCATTERING_MU_S * SCATTERING_NU * 3);

	//#pragma omp parallel for collapse(4)
	//for(int l = 0; l < SCATTERING_NU; l++)
	//	for(int k = 0; k < SCATTERING_MU_S; k++)
	//		for(int j = 0; j < SCATTERING_MU; j++)
	//			for (int i = 0; i < SCATTERING_R; i++) {
	//				float r = atmosphere.earthRad + (atmosphere.atmosphereRad - atmosphere.earthRad) * float(i) / (SCATTERING_R - 1);
	//				float mu = -1.0f + 2.0f * float(j) / (SCATTERING_MU - 1);
	//				float muS = -1.0f + 2.0f * float(k) / (SCATTERING_MU_S - 1);
	//				float nu = -1.0f + 2.0f * float(l) / (SCATTERING_NU - 1);
	//				glm::vec3 scattering = computeScattering(r, mu, muS, nu);

	//				int idx = ((l * SCATTERING_MU_S * SCATTERING_MU * SCATTERING_R + k * SCATTERING_MU * SCATTERING_R + j * SCATTERING_R + i) * 3);
	//				
	//				data[idx + 0] = scattering.x;
	//				data[idx + 1] = scattering.y;
	//				data[idx + 2] = scattering.z;

	//			}

	for (u8 i = 1; i < 4; i++)
		CO[i].setup(SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT, SCATTERING_TEXTURE_DEPTH, i); // 1 = scatteringLUT, 2 = rayleighDeltaLUT, 3 = mieDeltaLUT
	CO[4].setup(IRRADIANCE_TEXTURE_WIDTH, IRRADIANCE_TEXTURE_HEIGHT, 4); // 4 = directIrradiance
	CO[5].setup(SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT, SCATTERING_TEXTURE_DEPTH, 5); // 5 = scatteringDesntiyLUT
	CO[6].setup(IRRADIANCE_TEXTURE_WIDTH, IRRADIANCE_TEXTURE_HEIGHT, 6); // 6 = indirectIrradiance
	CO[7].setup(SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT, SCATTERING_TEXTURE_DEPTH, 7); // 7 = multiScatteringLUT

}
//
//void atmosphereLUTs::uvtoTransmittanceParams(float u, float v, float& r, float& mu) {
//	float H = sqrt(atmosphere.atmosphereRad * atmosphere.atmosphereRad - atmosphere.earthRad * atmosphere.earthRad);
//	float rho = v * H;
//	r = sqrt(rho * rho + atmosphere.earthRad * atmosphere.earthRad);
//	float d_min = r - atmosphere.earthRad;
//	float d_max = rho + H;
//	float d = d_min + u * (d_max - d_min);
//	mu = d == 0.0f ? 1.0f : (H * H - rho * rho - d * d) / (2.0f * r * d);
//	mu = clamp(mu, -1.0f, 1.0f);
//}
//
//glm::vec3 atmosphereLUTs::computeTransmittance(float r, float mu) {
//	const int STEPS = 250;
//	float dt = rayIntersectSphere(r, mu, atmosphere.atmosphereRad);
//	dt /= float(STEPS);
//
//	glm::vec3 transmittance = glm::vec3(1.0f, 1.0f, 1.0f);
//	float t = 0.0f;
//
//	for (int i = 0; i < STEPS; i++) {
//		float r_sample = sqrt(r * r + t * t + 2.0f * r * mu * t);
//		float height = r_sample - atmosphere.earthRad;
//	
//		glm::vec3 extinction = atmosphere.betaR * exp(-height / atmosphere.Hr) + atmosphere.betaM * exp(-height / atmosphere.Hm);
//		transmittance *= exp(-extinction * dt);
//		t += dt;
//
//	}
//	return transmittance;
//}
//
//glm::vec3 atmosphereLUTs::computeScattering(float r, float mu, float mu_s, float nu) {
//	const int STEPS = 50;
//	float dt = rayIntersectSphere(r, mu, atmosphere.atmosphereRad);
//	dt /= float(STEPS);
//
//	glm::vec3 result = { 0.0f, 0.0f, 0.0f };
//	float t = 0.0f;
//
//	for(int i = 0; i < STEPS; i++) {
//		float r_sample = sqrt(r * r + t * t + 2.0 * r * mu * t);
//		glm::vec3 pos = { 0.0f, r_sample, 0.0f };
//		float height = r_sample - atmosphere.earthRad;
//
//		glm::vec3 rayleigh = atmosphere.betaR * exp(-height/ atmosphere.Hr);
//		glm::vec3 mie = atmosphere.betaM * exp(-height/ atmosphere.Hm);
//	
//		float cos_theta = mu * mu_s + sqrt((1.0f - mu * mu) * (1.0f - mu_s * mu_s)) * nu;
//		float rayleigh_phase = 3.0f / (16.0f * M_PI) * (1.0f + cos_theta * cos_theta);
//
//		float g = atmosphere.meiG;
//
//		float mie_phase = 3.0f / (8.0f * M_PI) * ((1.0f - g * g) * (1.0f + cos_theta * cos_theta)) / ((2.0f + g * g) * pow(1.0f + g * g - 2.0f * g * cos_theta, 1.5f));
//
//		glm::vec3 sun_transmittance = computeTransmittance(r_sample, mu_s);
//		glm::vec3 view_transmittance = computeTransmittance(r_sample, -mu);
//
//		result += (rayleigh * rayleigh_phase + mie * mie_phase) * sun_transmittance * view_transmittance * dt;
//
//		t += dt;
//	}
//	return result;
//}
//
//float atmosphereLUTs::rayIntersectSphere(float r, float mu, float radius) {
//	float b = 2.0f * r * mu;
//	float c = r * r - radius * radius;
//
//	float det = b * b - 4.0f * c;
//
//	if (det < 0.0f) return 0.0f;
//	det = sqrt(det);
//
//	float t1 = (-b + det) * 0.5f;
//	float t2 = (-b - det) * 0.5f;
//
//	return (t1 >= 0.0f) ? t1 : t2;
//}