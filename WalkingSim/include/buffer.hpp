#pragma once
#include "utils.hpp"

class buffer { //I WANT TO DELETE BUFFERS AND SHIT FOR BETTER OPTIMSIATION
private:
	GLuint VAO;
	GLuint VBO[bufferID::NumBuffers];
	size_t num; 

	

public:
	buffer(std::vector<Vertex> vertices, drawFreq usage);
	buffer(std::vector<Vertex> vertices, std::vector<GLuint> indices, drawFreq usage);

	void draw(drawID dI, drawType dT, glm::vec3 pos,glm::vec3 rotation, glm::vec3 size, GLuint shaderID); //I'm expecting another funciton to call this after using textures
};