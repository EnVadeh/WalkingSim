#pragma once
#include "utils.hpp"

class textureManager {
private:
	std::vector<GLuint> textures;
	size_t texNum = 0;
	std::vector<std::string> texNames;
public:
	void loadTexture(const std::string& path, const std::string& name);//WILL OVERLOAD FOR CUBEMAP
	void bindTexture(size_t unit, size_t index,  size_t count, GLuint shaderID);
};

class computeOutput {
private:
	GLuint texture;
	std::string name;
public:
	void setup(size_t width, size_t height, u8 unit);
	void setup(size_t width, size_t height, size_t depth, u8 unit);
	void bind(GLuint shaderID, u8 unit, std::string name);
};