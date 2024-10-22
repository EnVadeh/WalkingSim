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
	glBindProgramPipeline(shaderID);
	glm::mat4 matModel = createGeometricToWorldMatrix(pos, rotation, size);
	setUniform(shaderID, "matModel", matModel);
	glPatchParameteri(GL_PATCH_VERTICES, 4);
	switch(dI) {
		case drawID::arrayDraw:
			glDrawArrays(dT, 0, num);
			break;
		case drawID::elementDraw:
			glDrawElements(dT, num, GL_UNSIGNED_INT, 0);
	}
}