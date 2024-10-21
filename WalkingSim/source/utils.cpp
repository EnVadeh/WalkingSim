#include "utils.hpp"

size_t SHADER_COUNT = 0;
std::vector <GLuint> SHADERS;
float deltaTime = 0.0f; // Time between current frame and last frame
float lastFrame = 0.0f;

glm::mat4 createGeometricToWorldMatrix(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale) {
	glm::mat4 Scale = glm::scale(glm::mat4(1.0f), scale);
	glm::mat4 Trans = glm::translate(glm::mat4(1.0f), position);
	return Trans * Scale;
}

//don't know weather to inline these funciton sin the hpp or just implement ithere
void setUniform(GLuint shaderID, std::string name, unsigned int val)
{
	const char* uniformNameptr = &name[0];
	GLint Loc = glGetUniformLocation(shaderID, uniformNameptr);
	glUniform1ui(Loc, val);
}
void setUniform(GLuint shaderID, std::string name, int val)
{
	const char* uniformNameptr = &name[0];
	GLint Loc = glGetUniformLocation(shaderID, uniformNameptr);
	glUniform1i(Loc, val);
}
void setUniform(GLuint shaderID, std::string name, float val)
{
	const char* uniformNameptr = &name[0];
	GLint Loc = glGetUniformLocation(shaderID, uniformNameptr);
	glUniform1f(Loc, val);
}
void setUniform(GLuint shaderID, std::string name, glm::vec3 val)
{
	const char* uniformNameptr = &name[0];
	GLint Loc = glGetUniformLocation(shaderID, uniformNameptr);
	glUniform3fv(Loc, 1, glm::value_ptr(val));
}
void setUniform(GLuint shaderID, std::string name, glm::vec4 val)
{
	const char* uniformNameptr = &name[0];
	GLint Loc = glGetUniformLocation(shaderID, uniformNameptr);
	glUniform4fv(Loc, 1, glm::value_ptr(val));
}
void setUniform(GLuint shaderID, std::string name, glm::mat3 val)
{
	const char* uniformNameptr = &name[0];
	GLint Loc = glGetUniformLocation(shaderID, uniformNameptr);
	glUniformMatrix3fv(Loc, 1, GL_FALSE, glm::value_ptr(val));
}
void setUniform(GLuint shaderID, std::string name, glm::mat4 val) {
	const char* uniformNameptr = &name[0];
	GLint Loc = glGetUniformLocation(shaderID, uniformNameptr);
	glUniformMatrix4fv(Loc, 1, GL_FALSE, glm::value_ptr(val));
}