#include "texture.hpp"
#include "stb_image.h"
#include <omp.h>

void textureManager::loadTexture(const std::string& path, const std::string& name) {
	GLuint texID;
	glCreateTextures(GL_TEXTURE_2D, 1, &texID);

	glTextureParameteri(texID, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTextureParameteri(texID, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTextureParameteri(texID, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTextureParameteri(texID, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	
	int width, height, nrChannels;
	unsigned char* data = stbi_load(path.c_str(), &width, &height, &nrChannels, 0);
	
	if (data) {
		GLenum baseFormat, sizedFormat;
		if (nrChannels == 3) {
			baseFormat = GL_RGB;
			sizedFormat = GL_RGB8;
		}
		else {//IM ASSUMING ILL NEVER USE AN IMAGE WITH ONLY ONE NR CHANNEL
			baseFormat = GL_RGBA;
			sizedFormat = GL_RGBA8;
		}
		glTextureStorage2D(texID, 1, sizedFormat, width, height);
		glTextureSubImage2D(texID, 0, 0, 0, width, height, baseFormat, GL_UNSIGNED_BYTE, data);
		glGenerateTextureMipmap(texID);
		stbi_image_free(data);
		textures.push_back(texID);
		texNames.push_back(name);
		texNum++;
	}
	else
		std::cout << "Failed to load texture: " << path << std::endl;
}

void textureManager::bindTexture(size_t unit, size_t index, size_t count, GLuint shaderID) {
	if (index > texNum || index + count > texNum) {
		std::cout << "Texture out of range!" << std::endl;
		return;
	}
	for (size_t i = 0; i < count; i++) {
		glBindTextureUnit(unit + i, textures[index + i]);
		setUniform(shaderID, texNames[index + i], static_cast<int>(unit + i));
	}
}

void computeOutput::setup(size_t width, size_t height) {
	glCreateTextures(GL_TEXTURE_2D, 1, &texture);
	glTextureParameteri(texture, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);//don't want lookup table values to repeat!
	glTextureParameteri(texture, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTextureParameteri(texture, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTextureParameteri(texture, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTextureStorage2D(texture, 1, GL_RGBA16F, width, height);
	glBindImageTexture(0, texture, 0, GL_FALSE, 0, GL_READ_WRITE, GL_RGBA16F);
}

void computeOutput::bind(GLuint shaderID) {
	glBindTextureUnit(0, texture);
	int x = 0;
	setUniform(shaderID, "LUT", x);
}