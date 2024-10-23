#pragma once
#include "utils.hpp"
#include "buffer.hpp"

struct light {
	glm::vec3 vLightPos = glm::vec3(0, 0, 0);
	glm::vec3 vLightDir = glm::vec3(0, 0, 0);
	glm::vec3 vLightColor = glm::vec3(0, 0, 0);
	glm::vec3 vEndPointA = glm::vec3(0, 0, 0);
	glm::vec3 vEndPointB = glm::vec3(0, 0, 0);
	bool isOn;
	lightType type;
};


struct lightUBOData {
	glm::vec3 vLightPos = glm::vec3(0, 0, 0);
	glm::vec3 vLightDir = glm::vec3(0, 0, 0);
	glm::vec3 vLightColor = glm::vec3(0, 0, 0);
	glm::vec3 vEndPointA = glm::vec3(0, 0, 0);
	glm::vec3 vEndPointB = glm::vec3(0, 0, 0);
};

class lightManager {
private:
	std::vector<light> lights;
	lightType eType;
	std::vector<lightUBOData> lightData;
public:
	void initLight(glm::vec3 lightPos, glm::vec3 lightDir, glm::vec3 lightColor);
	void initLight(glm::vec3 lightDir, glm::vec3 lightColor);
	void initLight(glm::vec3 lightPos, glm::vec3 lightDir, glm::vec3 lightColor, glm::vec3 epA, glm::vec3 epB);
	void turnOn(unsigned int lightIndex);
	void turnOn(unsigned int from, unsigned int to);
	void turnOff(unsigned int lightIndex);
	void turnOff(unsigned int from, unsigned int to);
	void setLights();
};