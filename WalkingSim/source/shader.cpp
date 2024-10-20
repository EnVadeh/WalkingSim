#include "shader.hpp"

shader::shader(std::string filepath1, std::string filepath2) {
	tess = false;

	std::ifstream shaderFile1(filepath1);
	std::string line1;
	while (std::getline(shaderFile1, line1)) {
		simpleShader.vertexSource += line1 + '\n';
	}
	shaderFile1.close();

	std::ifstream shaderFile2(filepath2);
	std::string line2;
	while (std::getline(shaderFile2, line2)) {
		simpleShader.fragmentSource += line2 + '\n';
	}
	shaderFile2.close();

}

shader::shader(std::string filepath1, std::string filepath2, std::string filepath3, std::string filepath4) {
	tess = true;

	std::ifstream shaderFile1(filepath1);
	std::string line1;
	while (std::getline(shaderFile1, line1)) {
		tessShader.vertexSource += line1 + '\n';
	}
	shaderFile1.close();

	std::ifstream shaderFile2(filepath2);
	std::string line2;
	while (std::getline(shaderFile2, line2)) {
		tessShader.tcSource += line2 + '\n';
	}
	shaderFile2.close();

	std::ifstream shaderFile3(filepath3);
	std::string line3;
	while (std::getline(shaderFile3, line3)) {
		tessShader.teSource += line3 + '\n';
	}
	shaderFile3.close();

	std::ifstream shaderFile4(filepath4);
	std::string line4;
	while (std::getline(shaderFile4, line4)) {
		tessShader.fragmentSource += line4 + '\n';
	}
	shaderFile4.close();

}

GLuint shader::compileShader(GLuint shaderType, const std::string& shaderSource) {
	GLuint shaderID = glCreateShader(shaderType);
	const char* src = shaderSource.c_str();
	glShaderSource(shaderID, 1, &src, NULL); //just like vaos and vbos, how many and what offset? 1 and NULL
	glCompileShader(shaderID);
	return shaderID;
		
}

GLuint shader::createShader() {
	program = glCreateProgram();
	if (!tess) {
		GLuint vs = compileShader(GL_VERTEX_SHADER, simpleShader.vertexSource);
		GLuint fs = compileShader(GL_FRAGMENT_SHADER, simpleShader.fragmentSource);
		glAttachShader(program, vs);
		glAttachShader(program, fs);
	}
	else {
		GLuint vs = compileShader(GL_VERTEX_SHADER, tessShader.vertexSource);
		GLuint tcs = compileShader(GL_TESS_CONTROL_SHADER, tessShader.tcSource);
		GLuint tes = compileShader(GL_TESS_EVALUATION_SHADER, tessShader.teSource);
		GLuint fs = compileShader(GL_FRAGMENT_SHADER, tessShader.fragmentSource);
		glAttachShader(program, vs);
		glAttachShader(program, tcs);
		glAttachShader(program, tes);
		glAttachShader(program, fs);
	}

	glLinkProgram(program);
	glValidateProgram(program);
	GLint infoLogLength = 512;
	int success;
	char infoLog[512];
	glGetProgramiv(program, GL_LINK_STATUS, &success);
	if (!success)
	{
		glGetProgramInfoLog(program, 512, NULL, infoLog);
		std::cout << "ERROR::SHADER::PROGRAM::LINKING_FAILED\n" << infoLog << std::endl;
	}

	return program;
}

