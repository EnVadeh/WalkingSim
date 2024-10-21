#pragma once
#include <vector>
#include <string>
#include <GL/glew.h>
#include <glm/exponential.hpp>
#include "stb_image.h"
#include <glm/gtc/matrix_transform.hpp>
#include <glm/ext.hpp>

struct bufferAttribs {
	//enum VAO_IDs { Triangles, NumVAOs };
	enum bufferIDs { ArrayBuffer, ElementBuffer, NumBuffers };
	enum attribIDs { vPos, vTex, vNormal, vTangent, vBiTangent };
	enum ssBuffers { vAlbedo, NumSSBs };
	enum ssBufferSizes { VertexSSBO = 16000000, floatSSBO = 4000000, vec3SSBO = 12000000 }; //i shoudl probably get size for cstom like if 1000 vertices: 1000 * sizeof(glm::vec4) type shit
	//GLuint Buffers[NumBuffers];
	//GLuint VAOs[NumVAOs];
};

static const std::string renderTextures[3] = {
	"colorRT",
	"normalRT",
	"depthRT"
};

static size_t SHADER_COUNT = 0;

std::vector <GLuint> SHADERS;

void setUniform(GLuint shaderID, std::string name, unsigned int val);
void setUniform(GLuint shaderID, std::string name, int val);
void setUniform(GLuint shaderID, std::string name, float val);
void setUniform(GLuint shaderID, std::string name, glm::vec3 val);
void setUniform(GLuint shaderID, std::string name, glm::vec4 val);
void setUniform(GLuint shaderID, std::string name, glm::mat3 val);
void setUniform(GLuint shaderID, std::string name, glm::mat4 val);