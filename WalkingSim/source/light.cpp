#include "light.hpp"

void lightManager::initLight(glm::vec3 lightPos, glm::vec3 lightDir, glm::vec3 lightColor) {
	light temp;
	temp.vLightPos = lightPos;
	temp.vLightDir = lightDir;
	temp.vLightColor = lightColor;
	temp.isOn = false;
	temp.type = lightType::pointLight;
	lights.push_back(temp);
	LIGHT_COUNT++;
}

void lightManager::initLight(glm::vec3 lightDir, glm::vec3 lightColor) {
	light temp;
	temp.vLightDir = lightDir;
	temp.vLightColor = lightColor;
	temp.isOn = false;
	temp.type = lightType::dirLight;
	lights.push_back(temp);
	LIGHT_COUNT++;
}

void lightManager::initLight(glm::vec3 lightPos, glm::vec3 lightDir, glm::vec3 lightColor, glm::vec3 epA, glm::vec3 epB) {
	light temp;
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
			lightData.push_back({ light.vLightPos, light.vLightDir, light.vLightColor, light.vEndPointA, light.vEndPointB });

	uniformBuffer<lightUBOData> lightBuff(lightData, drawFreq::staticDraw);
	lightBuff.bind();
}

