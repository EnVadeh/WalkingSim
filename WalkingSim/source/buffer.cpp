#include "buffer.hpp"

buffer::buffer(std::vector<Vertex> vertices, drawFreq usage) {
	glCreateVertexArrays(1, &VAO);
	glCreateBuffers(bufferID::NumBuffers, VBO);
	glNamedBufferData(VBO[bufferID::ArrayBuffer], sizeof(Vertex) * vertices.size(), vertices.data(), usage);

	glVertexArrayAttribFormat(VAO, attribID::vPos, 3, GL_FLOAT, GL_FALSE, offsetof(Vertex, vPosition));
	glEnableVertexArrayAttrib(VAO, attribID::vPos);
	glVertexArrayVertexBuffer(VAO, attribID::vPos, VBO[bufferID::ArrayBuffer], 0, sizeof(Vertex));

	glVertexArrayAttribFormat(VAO, attribID::vTex, 2, GL_FLOAT, GL_FALSE, offsetof(Vertex, vTex));
	glEnableVertexArrayAttrib(VAO, attribID::vTex);
	glVertexArrayVertexBuffer(VAO, attribID::vTex, VBO[bufferID::ArrayBuffer], 0, sizeof(Vertex));

	glVertexArrayAttribFormat(VAO, attribID::vNormal, 3, GL_FLOAT, GL_FALSE, offsetof(Vertex, vNormal));
	glEnableVertexArrayAttrib(VAO, attribID::vNormal);
	glVertexArrayVertexBuffer(VAO, attribID::vNormal, VBO[bufferID::ArrayBuffer], 0, sizeof(Vertex));

	num = vertices.size();
}


buffer::buffer(std::vector<Vertex> vertices, std::vector<GLuint> indices, drawFreq usage) {

	glCreateVertexArrays(1, &VAO);
	glCreateBuffers(bufferID::NumBuffers, VBO);
	glNamedBufferData(VBO[bufferID::ArrayBuffer], sizeof(Vertex) * vertices.size(), vertices.data(), usage);
	glNamedBufferData(VBO[bufferID::ElementBuffer], sizeof(indices), indices.data(), usage);

	glVertexArrayAttribFormat(VAO, attribID::vPos, 3, GL_FLOAT, GL_FALSE, offsetof(Vertex, vPosition));
	glEnableVertexArrayAttrib(VAO, attribID::vPos);
	glVertexArrayVertexBuffer(VAO, attribID::vPos, VBO[bufferID::ArrayBuffer], 0, sizeof(Vertex));

	glVertexArrayAttribFormat(VAO, attribID::vTex, 2, GL_FLOAT, GL_FALSE, offsetof(Vertex, vTex));
	glEnableVertexArrayAttrib(VAO, attribID::vTex);
	glVertexArrayVertexBuffer(VAO, attribID::vTex, VBO[bufferID::ArrayBuffer], 0, sizeof(Vertex));

	glVertexArrayAttribFormat(VAO, attribID::vNormal, 3, GL_FLOAT, GL_FALSE, offsetof(Vertex, vNormal));
	glEnableVertexArrayAttrib(VAO, attribID::vNormal);
	glVertexArrayVertexBuffer(VAO, attribID::vNormal, VBO[bufferID::ArrayBuffer], 0, sizeof(Vertex));

	glVertexArrayElementBuffer(VAO, VBO[bufferID::ElementBuffer]);
	num = indices.size();
}

void buffer::draw(drawID dI, drawType dT, glm::vec3 pos, glm::vec3 rotation,  glm::vec3 size, GLuint shaderID) {
	glBindVertexArray(VAO);
	glUseProgram(shaderID);
	glm::mat4 matModel = createGeometricToWorldMatrix(pos, rotation, size);
	setUniform(shaderID, "matModel", matModel);
	glPatchParameteri(GL_PATCH_VERTICES, 4);
	switch(dI) {
		case drawID::arrayDraw:
			glDrawArrays(dT, 0, num);
			break;
		case drawID::elementDraw:
			glDrawElements(dT, num, GL_UNSIGNED_INT, NULL);
			break;
	}
	
}

frameBuffer::frameBuffer() {
	glCreateFramebuffers(1, &FBO);
	std::cout << "The consturcotr is being called!";
	RT.resize(3);
	glCreateTextures(GL_TEXTURE_2D, 3, RT.data());
	for (size_t i = 0; i < 3; i++) {
		glTextureStorage2D(RT[i], 1, GL_RGB8, 2048, 2048);
		
		glTextureParameteri(RT[i], GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTextureParameteri(RT[i], GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTextureParameteri(RT[i], GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTextureParameteri(RT[i], GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	
		glNamedFramebufferTexture(FBO, GL_COLOR_ATTACHMENT0 + i, RT[i], 0);
	}
	setupRenderBuffer();
	glNamedFramebufferDrawBuffers(FBO, 3, drawBuffers);
}

frameBuffer::frameBuffer(bool depthOnly) : depthOnly(depthOnly) {
	glCreateFramebuffers(1, &FBO);

	glCreateTextures(GL_TEXTURE_2D, 1, &shadowRT);

	glTextureStorage2D(shadowRT, 1, GL_DEPTH_COMPONENT, 2048, 2048);
	
	glTextureParameteri(shadowRT, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTextureParameteri(shadowRT, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTextureParameteri(shadowRT, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTextureParameteri(shadowRT, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	glNamedFramebufferTexture(FBO, GL_DEPTH_ATTACHMENT, shadowRT, 0);

	glNamedFramebufferDrawBuffer(FBO, GL_NONE);
}

void frameBuffer::setupRenderBuffer() {
	glCreateRenderbuffers(1, &RBO);
	glNamedRenderbufferStorage(RBO, GL_DEPTH_STENCIL, 2048, 2048);
	glNamedFramebufferRenderbuffer(FBO, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, RBO);

}
void frameBuffer::bind() {
	glViewport(0, 0, 2048, 2048);
	glBindFramebuffer(GL_FRAMEBUFFER, FBO);
	if (!depthOnly)
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
	else
		glClear(GL_DEPTH_BUFFER_BIT);
}

void frameBuffer::sample(GLuint shaderID, GLint base_unit) {
	if (depthOnly) {
		glBindTextureUnit(base_unit, shadowRT);
		setUniform(shaderID, "shadowRT", base_unit);
	}
	else {
		std::vector<std::string> texNames;
		texNames.push_back("colorRT");
		texNames.push_back("normalRT");
		texNames.push_back("depthRT");
		for (int i = 0; i < 3; i++) {
			glBindTextureUnit(base_unit + i, RT[i]);
			setUniform(shaderID, texNames[i], base_unit + i);
		}
	}

}

terrain::terrain(size_t length, size_t breadth) : length(length), breadth(breadth) { 
	std::vector<Vertex> tVertices;
	for(size_t i = 0; i < length; ++i)
		for (size_t j = 0; j < breadth; ++j) {
			Vertex temp;
			temp.vPosition = glm::vec3(i, 0, j);
			temp.vNormal = glm::vec3(0, 1, 0); //cause it's flat
			temp.vTex = glm::vec2(float(i) / (length - 1), float(j) / (breadth - 1));
			tVertices.emplace_back(temp);
			//std::cout << "vertices generated are: (" << temp.vPosition.x << ", " << temp.vPosition.y << ", " << temp.vPosition.z << ")" << std::endl;
			//std::cout << "normals generated are: (" << temp.vNormal.x << ", " << temp.vNormal.y << ", " << temp.vNormal.z << ")" << std::endl;
			//std::cout << "texcoord generated are: (" << temp.vTex.x << ", " << temp.vTex.y <<")" << std::endl;
		}

	std::vector<GLuint> indices;
	for (int i = 0; i < length; ++i) {
		for (int j = 0; j < breadth; ++j) {
			int start = i * breadth + j;
			indices.push_back(start);
			indices.push_back(start + 1);
			indices.push_back(start + breadth);

			indices.push_back(start + 1);
			indices.push_back(start + breadth + 1);
			indices.push_back(start + breadth);
		}
	}
	tBuffer = std::make_shared<buffer>(tVertices, indices, drawFreq::staticDraw);
}

void terrain::draw(GLuint shaderID, glm::vec3 pos, glm::vec3 size) {

	tBuffer->draw(elementDraw, triDraw, pos, glm::vec3(0, 0, 0), size, shaderID);
}

screenQuad::screenQuad(size_t length, size_t breadth) : length(length), breadth(breadth) {
	std::vector<Vertex> quadVertices = {
		// Positions					// Texture Coords		// Normals (0, 0, 0)
		{ glm::vec3(-1.0f,  1.0f, 0.0f), glm::vec2(0.0f, 1.0f), glm::vec3(0.0f, 0.0f, 0.0f) },
		{ glm::vec3(-1.0f, -1.0f, 0.0f), glm::vec2(0.0f, 0.0f), glm::vec3(0.0f, 0.0f, 0.0f) },
		{ glm::vec3(1.0f, -1.0f, 0.0f), glm::vec2(1.0f, 0.0f), glm::vec3(0.0f, 0.0f, 0.0f) },

		{ glm::vec3(-1.0f,  1.0f, 0.0f), glm::vec2(0.0f, 1.0f), glm::vec3(0.0f, 0.0f, 0.0f) },
		{ glm::vec3(1.0f, -1.0f, 0.0f), glm::vec2(1.0f, 0.0f), glm::vec3(0.0f, 0.0f, 0.0f) },
		{ glm::vec3(1.0f,  1.0f, 0.0f), glm::vec2(1.0f, 1.0f), glm::vec3(0.0f, 0.0f, 0.0f) }
	};

	qBuffer = std::make_shared<buffer>(quadVertices, drawFreq::staticDraw);
}

void screenQuad::draw(GLuint shaderID) {
	glViewport(0, 0, length, breadth);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
	qBuffer->draw(arrayDraw, triDraw, glm::vec3(0, 0, 0), glm::vec3(0, 0, 0), glm::vec3(1, 1, 1), shaderID);
}