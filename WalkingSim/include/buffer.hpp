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

class terrain{
private: //use these while positioning to make sure that the position is always "in the middle"
	size_t length; //s
	size_t breadth; //t
	std::shared_ptr<buffer> tBuffer;

public:
	terrain(size_t length, size_t breadth);
	void draw(GLuint shaderID, glm::vec3 pos, glm::vec3 size);  //I wanna take in camera position so I can load chunks later I think.. But not in this file

};