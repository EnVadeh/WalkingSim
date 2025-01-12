#pragma once
#include <vector>
#include <string>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/exponential.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/ext.hpp>
#include <iostream>
#include <memory>
#include <chrono>
#include "camera.hpp"
#include <cstdint>

#define M_PI 3.14159265358979323846
using u8 = uint8_t;

//enum VAO_IDs { Triangles, NumVAOs };
enum bufferID { ArrayBuffer, ElementBuffer, NumBuffers };
enum attribID { vPos, vTex, vNormal, vTangent, vBiTangent };
enum ssBuffer { vAlbedo, NumSSBs };
enum ssBufferSize { VertexSSBO = 16000000, floatSSBO = 4000000, vec3SSBO = 12000000 }; //i shoudl probably get size for cstom like if 1000 vertices: 1000 * sizeof(glm::vec4) type shit
enum drawID { arrayDraw, elementDraw };
enum drawFreq { staticDraw = GL_STATIC_DRAW, dynamicDraw = GL_DYNAMIC_DRAW };
enum drawType { triDraw = GL_TRIANGLES, patchDraw = GL_PATCHES };
enum lightType { dirLight, pointLight, areaLight };


struct Vertex { //THIS IS HOW I WANT EVERYTHING TO BE DEFINED!
	glm::vec3 vPosition;
	glm::vec2 vTex;
	glm::vec3 vNormal;
};


static const std::string renderTextures[3] = {
	"colorRT",
	"normalRT",
	"depthRT"
};

struct densityProfileLayer {
	float width;
	float exp_term;
	float exp_scale;
	float linear_term;
	float constant_term;
	float padding1;
	float padding2;
	float padding3;
};


template <typename T>
T clamp(T value, T minVal, T maxVal) {
    if (value < minVal) {
        return minVal;
    }
    else if (value > maxVal) {
        return maxVal;
    }
    else {
        return value;
    }
}


extern float deltaTime; // Time between current frame and last frame
extern float lastFrame;

extern size_t SHADER_COUNT;
extern size_t LIGHT_COUNT;
extern size_t UBO_COUNT;

extern std::vector <GLuint> SHADERS;

glm::mat4 createGeometricToWorldMatrix(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale);

void setCommonShader(const std::string& path);

void setUniform(GLuint shaderID, std::string name, unsigned int val);
void setUniform(GLuint shaderID, std::string name, int val);
void setUniform(GLuint shaderID, std::string name, float val);
void setUniform(GLuint shaderID, std::string name, glm::vec3 val);
void setUniform(GLuint shaderID, std::string name, glm::vec4 val);
void setUniform(GLuint shaderID, std::string name, glm::mat3 val);
void setUniform(GLuint shaderID, std::string name, glm::mat4 val);

template <typename T>
class image2D{ //I don't think it'd matter if UV starts from the 0 of the vector or the (width * height) - width, because I'm the one putting values anyways
private:
	std::vector<T> image;
	std::vector<bool> cFlag;
	size_t width;
	size_t height;
public:
	image2D(size_t row, size_t columns);
	void write(T data, size_t x, size_t y, bool replace);
	T read(size_t x, size_t y);
	glm::vec2 size();
};

template <typename T>
image2D<T>::image2D(size_t row, size_t columns) {
	width = row;
	height = columns;
	image.resize(width * height);
	cFlag.resize(width * height, false);
}

template <typename T>
void image2D<T>::write(T data, size_t x, size_t y, bool replace){
	if (x >= width || y >= height) {
		std::cout << "pixel out of bounds, the size of the image is: " << width << " , " << height << std::endl;
		return;
	}
	if (!replace)
		if (cFlag[x + y * width] == false) {
			image[x + y * width] = data;
			cFlag[x + y * width] = true;
			return;
		}
		else {
			std::cout << "Already written" << std::endl; //remember that the width goes from 0 -> width-1
			return;
		}
	image[x + y * width] = data; 
	cFlag[x + y * width] = true;
	return;
}

template <typename T>
T image2D<T>::read(size_t x, size_t y)
{
	if (x >= width || y >= height) {
		std::cout << "pixel out of bounds, the size of the image is: " << width << " , " << height << std::endl;
		return static_cast<T>(0);
	}
	if (cFlag[x + y * width] == false) {
		return static_cast<T>(0);
	}
	return image[x + y * width];
}

template <typename T>
glm::vec2 image2D<T>::size() { return glm::vec2{ width, height }; }