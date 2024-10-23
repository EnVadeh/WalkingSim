#pragma once
#include "utils.hpp"
#include "buffer.hpp"

struct light {
	glm::vec4 vLightPos = glm::vec4(0, 0, 0, 0);
	glm::vec4 vLightDir = glm::vec4(0, 0, 0, 0);
	glm::vec4 vLightColor = glm::vec4(0, 0, 0, 0);
	glm::vec4 vEndPointA = glm::vec4(0, 0, 0, 0);
	glm::vec4 vEndPointB = glm::vec4(0, 0, 0, 0);
	bool isOn;
	lightType type;
};


struct lightUBOData {
	glm::vec4 vLightPos = glm::vec4(0, 0, 0, 0);
	glm::vec4 vLightDir = glm::vec4(0, 0, 0, 0);
	glm::vec4 vLightColor = glm::vec4(0, 0, 0, 0);
	glm::vec4 vEndPointA = glm::vec4(0, 0, 0, 0);
	glm::vec4 vEndPointB = glm::vec4(0, 0, 0, 0);
};

class lightManager {
private:
	std::vector<light> lights;
	lightType eType;
	std::vector<lightUBOData> lightData;
public:
	void initLight(glm::vec4 lightPos, glm::vec4 lightDir, glm::vec4 lightColor);
	void initLight(glm::vec4 lightDir, glm::vec4 lightColor);
	void initLight(glm::vec4 lightPos, glm::vec4 lightDir, glm::vec4 lightColor, glm::vec4 epA, glm::vec4 epB);
	void turnOn(unsigned int lightIndex);
	void turnOn(unsigned int from, unsigned int to);
	void turnOff(unsigned int lightIndex);
	void turnOff(unsigned int from, unsigned int to);
	void setLights();
};