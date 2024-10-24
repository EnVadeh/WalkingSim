#pragma once
#include "utils.hpp"

//do I want to passw by reference? and template the buffer too? let's see

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


class frameBuffer { //Should I create a vector of all RenderTargets? or like a struct of 3 rendertargets along with the FrameBuffer number, like which framebuffer invoked type shit
private:
	GLenum drawBuffers[3] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2 };
	GLuint FBO;
	GLuint RBO;
	std::vector<GLuint> RT;
	GLuint shadowRT;
	bool depthOnly = false;
	void setupRenderBuffer(); //I shall help make different types of frameBuffers that use RBO vs RTs
public:
	frameBuffer(); 
	frameBuffer(bool depthOnly);
	void bind();
	void sample(GLuint shaderID, GLint base_unit); //the textures will be sent to shaders
};

template <typename T>
class uniformBuffer {
private:
	GLuint UBO;
	GLuint bindingPoint;
public:
	uniformBuffer(const std::vector<T>& data, drawFreq usage);
	int returnBP();
	void bind();
};

template <typename T>
uniformBuffer<T>::uniformBuffer(const std::vector<T>& data, drawFreq usage) {
	glCreateBuffers(1, &UBO);
	glNamedBufferData(UBO, data.size() * sizeof(T), data.data(), usage);
	bindingPoint = UBO_COUNT;
	UBO_COUNT++;
}

template<typename T>
int uniformBuffer<T>::returnBP() {
	return bindingPoint;
}

template <typename T>
void uniformBuffer<T>::bind() {
	glBindBufferBase(GL_UNIFORM_BUFFER, bindingPoint, UBO);
}

class terrain{
private: //use these while positioning to make sure that the position is always "in the middle"
	size_t length; //s
	size_t breadth; //t
	std::shared_ptr<buffer> tBuffer;

public:
	terrain(size_t length, size_t breadth);
	void draw(GLuint shaderID, glm::vec3 pos, glm::vec3 size);  //I wanna take in camera position so I can load chunks later I think.. But not in this file
};

class screenQuad {
private:
	std::shared_ptr<buffer> qBuffer;
	size_t length;
	size_t breadth;
public:
	screenQuad(size_t length, size_t breadth);
	void draw(GLuint shaderID);
};