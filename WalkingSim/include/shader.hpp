#pragma once
#include "utils.hpp"
#include <string>
#include <iostream>
#include <fstream>

struct computeShaderSource {
	std::string shaderSource;
};

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

class shader { //SHADSERS WONT BE AFFECTED MUCH BY DSA OR NOT
private:
	static std::string readFile(const std::string& filepath);
	computeShaderSource computeShader;
	simpleShaderSource simpleShader; //need better name
	tessShaderSource tessShader;
	bool tess;
	bool compute = 0;
	GLuint program;
	GLuint compileShader(GLuint shaderType, const std::string& shaderSource);
public:
	shader(const std::string& filepath);
	shader(const std::string& filepath1, const std::string& filepath2);//AM PASSING REFERENCES CHECK AGAIN
	shader(const std::string& filepath1, const std::string& filepath2, const std::string& filepath3, const std::string& filepath4);
	GLuint createShader(); //REMEMBER THIS WILL ADD ALL SHADERS TO AN ARRAY SO THAT THE CAMERA MATRICES ARE SENT
};

