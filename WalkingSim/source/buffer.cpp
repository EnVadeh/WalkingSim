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

