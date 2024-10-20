#pragma once
#include "utils.hpp"
#include <string>
#include <iostream>
#include <fstream>

struct simpleShaderSource {
	std::string vertexSource;
	std::string fragmentSource;
};

struct tessShaderSource {
	std::string vertexSource;
	std::string tcSource;
	std::string teSource;
	std::string fragmentSource;
};

class shader {
private:
	simpleShaderSource simpleShader; //need better name
	tessShaderSource tessShader;
	bool tess;
	GLuint program;
	GLuint compileShader(GLuint shaderType, const std::string& shaderSource);
public:
	shader(std::string filepath1, std::string filepath2);
	shader(std::string filepath1, std::string filepath2, std::string filepath3, std::string filepath4);
	GLuint createShader();
};