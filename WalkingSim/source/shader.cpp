#include "shader.hpp"

std::string shader::readFile(const std::string& filepath){
	std::ifstream file(filepath);
	if (!file.is_open()) {
		throw std::runtime_error("Failed to open shader file: " + filepath);
	}

	return std::string(
		std::istreambuf_iterator<char>(file),  //points to start of string buffer
		std::istreambuf_iterator<char>()	   //points to end of string buffer
	);
}

shader::shader(const std::string& filepath1, const std::string& filepath2) : tess(false) {

	simpleShader.vertexSource = readFile(filepath1);
	simpleShader.fragmentSource = readFile(filepath2);

}

shader::shader(const std::string& filepath1, const std::string& filepath2, const std::string& filepath3, const std::string& filepath4) : tess(true) {

	tessShader.vertexSource = readFile(filepath1);
	tessShader.tcSource = readFile(filepath2);
	tessShader.teSource = readFile(filepath3);
	tessShader.fragmentSource = readFile(filepath4);


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
		glLinkProgram(program);
		glValidateProgram(program);
		glDeleteShader(vs);
		glDeleteShader(fs);

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
		glLinkProgram(program);
		glValidateProgram(program);
		glDeleteShader(vs);
		glDeleteShader(tcs);
		glDeleteShader(tes);
		glDeleteShader(fs);
	}

	
	GLint infoLogLength = 512;
	int success;
	char infoLog[512];
	glGetProgramiv(program, GL_LINK_STATUS, &success);
	if (!success)
	{
		glGetProgramInfoLog(program, 512, NULL, infoLog);
		std::cout << "ERROR::SHADER::PROGRAM::LINKING_FAILED\n" << infoLog << std::endl;
	}
	SHADER_COUNT++;
	SHADERS.push_back(program);
	return program;
}