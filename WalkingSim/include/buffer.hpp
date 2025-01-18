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
	size_t breadth; //t
	std::shared_ptr<buffer> tBuffer;

public:
	size_t length; //s
	terrain(size_t length, size_t breadth);
	void draw(GLuint shaderID, glm::vec3 pos, glm::vec3 size);  //I wanna take in camera position so I can load chunks later I think.. But not in this file
};

class skyBox {
private:
	std::shared_ptr<buffer> sBuffer;
public:
	skyBox();
	void draw(GLuint shaderID);
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

//Let's try to use lambas I guess. So first we have a calculatePerlin function that takes in either Horizontal Lerp func or Vertical Lerp func? 

//This is broken, I need to make small grids here...
template<typename T>
void perlinNoise(image2D<T>& image) {
	//We are going to 'grid' the original image, where each grid is 50x50
	//We are going to store, only the 50x50 random vectors

	size_t size_x = image.size().x;
	size_t size_y = image.size().y;

	glm::vec2 cornerVecSize = { size_x / 100, size_y / 100};
	image2D<vector2D> cornerVecs(cornerVecSize.x, cornerVecSize.y);
	for (size_t x = 0; x < cornerVecSize.x; x++)
		for (size_t y = 0; y < cornerVecSize.y; y++)
			cornerVecs.write(x, y, randomGradient(vector2D(x, y)), 0);

	for (size_t x = 0; x < size_x; x++) {
		for (size_t y = 0; y < size_y; y++) {
			float u;
			float v;
			u = float(y % 100) / 100; //u goes towards 1 along the x-axis, and in the case of a 2D array, going across the x-axis is changing the value of y
			v = float(x % 100) / 100; //same as above but for v and y-axis

			size_t x_index = x / 100; //The indices of the abstract grid cells, which one we're in rn, it also starts from 0
			size_t y_index = y / 100;

			vector2D pos = vector2D(u, v);//current position in 0-1

			//remember that v goes towards 1 when going down...
			
			//direction vectors from the corners of a grid to the position
			vector2D tl = pos - vector2D(0,0);
			vector2D tr = pos - vector2D(1,0);
			vector2D bl = pos - vector2D(0,1);
			vector2D br = pos - vector2D(1,1);

			//Dot products between the corner vectors and the direction vectors from the corner to position
			double tld = dotProduct(tl, cornerVecs.read(x_index, y_index));
			double trd = dotProduct(tr, cornerVecs.read(x_index + 1, y_index));
			double bld = dotProduct(bl, cornerVecs.read(x_index, y_index + 1));
			double brd = dotProduct(br, cornerVecs.read(x_index + 1, y_index + 1));

			//interpolating the left and right vectors along top and bottom
			double t = mix(tld, trd, u);
			double b = mix(bld, brd, u);

			//interpolating the top and bottom vectors along the y axis
			double perlin = mix(b, t, v);
			
			//don't want y to be negative
			//perlin = fabs(perlin);
			//decreasing the value
			//perlin = perlin / 2;
			image.write(x, y, perlin, 0);
		}
	}
}