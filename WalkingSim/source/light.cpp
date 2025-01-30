#include "light.hpp"

void lightManager::initLight(glm::vec4 lightPos, glm::vec4 lightDir, glm::vec4 lightColor) {
	light temp;
	temp.vLightPos = lightPos;
	temp.vLightDir = lightDir;
	temp.vLightColor = lightColor;
	glm::vec3 vRight = glm::normalize(glm::cross(glm::vec3(lightDir), glm::vec3(0, 1, 0)));
	glm::vec3 vUp = glm::normalize(glm::cross(vRight, glm::vec3(lightDir)));
	temp.mView = glm::lookAt(glm::vec3(lightPos), glm::vec3(lightDir), vUp);
	temp.mProjection = glm::perspective(float(glm::radians(60.0f)), 1.0f, 0.1f, 10000.0f);
	temp.isOn = false;
	temp.type = lightType::pointLight;
	lights.push_back(temp);
	LIGHT_COUNT++;
}

void lightManager::initLight(glm::vec4 lightDir, glm::vec4 lightColor) { //I'd love to take a distance factor for lightPosition fake, and boundary values
	light temp;
	glm::vec3 lightPos = -glm::vec3(lightDir) * 5.0f;
	glm::vec3 vRight = glm::normalize(glm::cross(glm::vec3(lightDir), glm::vec3(0, 1, 0)));
	glm::vec3 vUp = glm::normalize(glm::cross(vRight, glm::vec3(lightDir)));
	temp.mView = glm::lookAt(glm::vec3(lightPos), glm::vec3(lightDir), glm::vec3(0, 1, 0));
	temp.mProjection = glm::ortho(-5.f, 5.f, -5.f, 5.f, -10000.1f, 10000.0f);
	temp.vLightDir = lightDir;
	temp.vLightColor = lightColor;
	temp.isOn = false;
	temp.type = lightType::dirLight;
	lights.push_back(temp);
	LIGHT_COUNT++;
}

void lightManager::initLight(glm::vec4 lightPos, glm::vec4 lightDir, glm::vec4 lightColor, glm::vec4 epA, glm::vec4 epB) {
	light temp;

	//I haven't figured out what to do with Area light Source
	//glm::vec3 vRight = glm::normalize(glm::cross(glm::vec3(lightDir), glm::vec3(0, 1, 0)));
	//glm::vec3 vUp = glm::normalize(glm::cross(vRight, glm::vec3(lightDir)));
	temp.vLightPos = lightPos;
	temp.vLightDir = lightDir;
	temp.vLightColor = lightColor;
	temp.isOn = false;
	temp.type = lightType::areaLight;
	temp.vEndPointA = epA;
	temp.vEndPointB = epB;
	lights.push_back(temp);
	LIGHT_COUNT++;
}

void lightManager::turnOn(unsigned int lightIndex) {
	if (lightIndex < LIGHT_COUNT)
		lights[lightIndex].isOn = true;
	else
		std::cout << "light " << lightIndex << " doesn't exist" << std::endl;
}

void lightManager::turnOn(unsigned int from, unsigned int to) {
	if (to < LIGHT_COUNT)
		for (unsigned int i = from; i < to; i++)
			lights[i].isOn = true;
	else
		std::cout << "only " << LIGHT_COUNT << " lights exist" << std::endl;
}

void lightManager::turnOff(unsigned int lightIndex) {
	if (lightIndex < LIGHT_COUNT)
		lights[lightIndex].isOn = false;
	else
		std::cout << "light " << lightIndex << " doesn't exist" << std::endl;
}

void lightManager::turnOff(unsigned int from, unsigned int to) {
	if (to < LIGHT_COUNT)
		for (unsigned int i = from; i < to; i++)
			lights[i].isOn = false;
	else
		std::cout << "only " << LIGHT_COUNT << " lights exist" << std::endl;
}

void lightManager::setLights() {
	for (const auto& light : lights)
		if (light.isOn)
			lightData.push_back({ light.vLightPos, light.vLightDir, light.vLightColor, light.vEndPointA, light.vEndPointB, light.mProjection, light.mView });

	uniformBuffer<lightUBOData> lightBuff(lightData, drawFreq::staticDraw);
	lightBuff.bind();
}

