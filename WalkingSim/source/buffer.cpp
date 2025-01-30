#include "buffer.hpp"
#include <fstream>

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
	glNamedBufferData(VBO[bufferID::ElementBuffer], sizeof(GLuint) * indices.size(), indices.data(), usage);

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
	RT.resize(3);
	glCreateTextures(GL_TEXTURE_2D, 3, RT.data());
	for (size_t i = 0; i < 3; i++) {
		glTextureStorage2D(RT[i], 1, GL_RGB16, 2048, 2048);
		
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
	RBO = -1;
	glCreateTextures(GL_TEXTURE_2D, 1, &shadowRT);

	glTextureStorage2D(shadowRT, 1, GL_DEPTH_COMPONENT32F, 4096, 4096);
	
	glTextureParameteri(shadowRT, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTextureParameteri(shadowRT, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTextureParameteri(shadowRT, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTextureParameteri(shadowRT, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	glNamedFramebufferTexture(FBO, GL_DEPTH_ATTACHMENT, shadowRT, 0);

	glNamedFramebufferDrawBuffer(FBO, GL_NONE);
	glNamedFramebufferReadBuffer(FBO, GL_NONE);
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
	//float stepSize = 0.03125;
	float stepSize = 0.015625; //runs about 1280 times me thinks
	float normalisingFactor = (float(length) + stepSize)/stepSize;
	size_t divisonNum = static_cast<size_t>(1.0f / stepSize);
	//generating HeightMap:

	size_t hm1Size = 500;
	size_t hm2Size = 200;
	image2D<double>heightMap1(hm1Size, hm1Size);
	image2D<double>heightMap2(hm2Size, hm2Size);
	perlinNoise(heightMap1, 50, true);
	perlinNoise(heightMap2, 20, false);
	float heightScale1 = 3.5;
	float heightScale2 = 0.8;

	//std::ofstream out("newfile.txt");
	//std::streambuf* coutbuf = std::cout.rdbuf(); // Save old buffer
	//std::cout.rdbuf(out.rdbuf());   // Redirect std::cout to file
	//std::cout << "P3\n";
	//std::cout << 500 << " " << 500 << "\n";
	//std::cout << "255\n";

	//for (size_t tempx = 0; tempx < 500; tempx++) {
	//	for (size_t tempy = 0; tempy < 500; tempy++) {
	//	//std::cout << static_cast<int>(heightMap.read(tempx, tempy) * 255) << " "<< static_cast<int>(heightMap.read(tempx, tempy) * 255) <<" "<< static_cast<int>(heightMap.read(tempx, tempy) * 255);
	//}
	//	//std::cout<<"\n";
	//}
	//Generating triangles row wise

	size_t perlinY = 0;
	for (float j = 0; j < breadth + stepSize; j = j + stepSize) {
	size_t perlinX = 0;
		for (float i = 0; i < length + stepSize; i = i + stepSize) {
			size_t hm1NormX = (float(perlinX)) / (normalisingFactor) * hm1Size;
			size_t hm1NormY = (float(perlinY)) / (normalisingFactor) * hm1Size;
			size_t hm2NormY = (float(perlinY)) / (normalisingFactor) * hm2Size;
			size_t hm2NormX = (float(perlinX)) / (normalisingFactor) * hm2Size;
			//I want to calculate the normal of this vertex, to do that I can basiaclly the cross product of vectors
			//from +x to -x and +z to -z
			//since the terrain is made from 0 to +x axis and 0 to + z axis:
			//less than 0 = fucked
			//glm::vec3 xVec = glm::vec3(perlinX + 3, heightMap1.clampRead(hm1NormX + 3, hm1NormY) + heightMap2.clampRead(hm2NormX + 3, hm2NormY), perlinY) - glm::vec3(perlinX > 3? perlinX - 3 : 0, heightMap1.clampRead(hm1NormX > 3? hm1NormX - 3 : 0, hm1NormY) + heightMap2.clampRead(hm2NormX > 3? hm2NormX - 3 : 0, hm2NormY) , perlinY);
			//glm::vec3 zVec = -glm::vec3(perlinX, heightMap1.clampRead(hm1NormX, hm1NormY + 3) + heightMap2.clampRead(hm2NormX, hm2NormY + 3), perlinY  + 3) + glm::vec3(perlinX, heightMap1.clampRead(hm1NormX, hm1NormY > 3 ? hm1NormY - 3 : 0) + heightMap2.clampRead(hm2NormX, hm2NormY > 3 ? hm2NormY - 3 : 0) , perlinY > 3? perlinY - 3: 0);

			float DX = (heightMap1.clampRead(hm1NormX + 3, hm1NormY) * heightScale1 + heightMap2.clampRead(hm2NormX + 3, hm2NormY) * heightScale2) - (heightMap1.clampRead(hm1NormX > 3 ? hm1NormX - 3 : 0, hm1NormY) * heightScale1 + heightMap2.clampRead(hm2NormX > 3 ? hm2NormX - 3 : 0, hm2NormY) * heightScale2);
			float DZ = -(heightMap1.clampRead(hm1NormX, hm1NormY > 3 ? hm1NormY - 3 : 0) * heightScale1 + heightMap2.clampRead(hm2NormX, hm2NormY > 3 ? hm2NormY - 3 : 0) * heightScale2) + (heightMap1.clampRead(hm1NormX, hm1NormY + 3) * heightScale1 + heightMap2.clampRead(hm2NormX, hm2NormY + 3) * heightScale2);
			glm::vec3 xVec = glm::vec3(1, DX, 0);
			glm::vec3 zVec = glm::vec3(0, DZ, 1);

			glm::vec3 norm = (glm::cross(zVec, xVec)) ;

			tVertices.push_back({ glm::vec3(0.0f + i, 0.0f + heightMap1.clampRead(hm1NormX, hm1NormY) * heightScale1 + heightMap2.safeRead(hm2NormX, hm2NormY) * heightScale2 , 0.0f + j), glm::vec2(i / length, j / breadth), norm});
			perlinX++;
		}
		perlinY++;
	}
	std::vector<GLuint> indices;

	//Generating indices quad wise which is row wise
	for (size_t j = 0; j < breadth * divisonNum ; j++) 
		for (size_t i = 0; i < length * divisonNum ; i++) {
			size_t i0 = i + j * (breadth * divisonNum + 1);	// top left corner
			size_t i1 = i + j * (breadth * divisonNum + 1) + 1;	// top right corner
			size_t i2 = i + (j + 1) * (breadth * divisonNum + 1);	// bottom left corner
			size_t i3 = i + (j + 1) * (breadth * divisonNum + 1) + 1;	// bottom right 
			//going to do the quad counter clockwise since the forward face is counterclock wise in my setup
			indices.push_back(i0);
			indices.push_back(i2);
			indices.push_back(i1);
			indices.push_back(i1);
			indices.push_back(i2);
			indices.push_back(i3);

	}

	tBuffer = std::make_shared<buffer>(tVertices, indices, drawFreq::staticDraw);
}

void terrain::draw(GLuint shaderID, glm::vec3 pos, glm::vec3 size) {
	tBuffer->draw(elementDraw, triDraw, pos, glm::vec3(0, 0, 0), size, shaderID);
}


skyBox::skyBox() {
	std::vector<Vertex> sVertices;

	sVertices = { 
	// Front face
	{{-1.0f,  1.0f, -1.0f}, {0.0f, 1.0f}, {0.0f, 0.0f, -1.0f}},
	{{-1.0f, -1.0f, -1.0f}, {0.0f, 0.0f}, {0.0f, 0.0f, -1.0f}},
	{{ 1.0f, -1.0f, -1.0f}, {1.0f, 0.0f}, {0.0f, 0.0f, -1.0f}},
	{{ 1.0f, -1.0f, -1.0f}, {1.0f, 0.0f}, {0.0f, 0.0f, -1.0f}},
	{{ 1.0f,  1.0f, -1.0f}, {1.0f, 1.0f}, {0.0f, 0.0f, -1.0f}},
	{{-1.0f,  1.0f, -1.0f}, {0.0f, 1.0f}, {0.0f, 0.0f, -1.0f}},

	// Back face
	{{-1.0f, -1.0f,  1.0f}, {0.0f, 0.0f}, {0.0f, 0.0f, 1.0f}},
	{{-1.0f,  1.0f,  1.0f}, {0.0f, 1.0f}, {0.0f, 0.0f, 1.0f}},
	{{ 1.0f,  1.0f,  1.0f}, {1.0f, 1.0f}, {0.0f, 0.0f, 1.0f}},
	{{ 1.0f,  1.0f,  1.0f}, {1.0f, 1.0f}, {0.0f, 0.0f, 1.0f}},
	{{ 1.0f, -1.0f,  1.0f}, {1.0f, 0.0f}, {0.0f, 0.0f, 1.0f}},
	{{-1.0f, -1.0f,  1.0f}, {0.0f, 0.0f}, {0.0f, 0.0f, 1.0f}},

	// Left face
	{{-1.0f,  1.0f,  1.0f}, {1.0f, 1.0f}, {-1.0f, 0.0f, 0.0f}},
	{{-1.0f, -1.0f,  1.0f}, {1.0f, 0.0f}, {-1.0f, 0.0f, 0.0f}},
	{{-1.0f, -1.0f, -1.0f}, {0.0f, 0.0f}, {-1.0f, 0.0f, 0.0f}},
	{{-1.0f, -1.0f, -1.0f}, {0.0f, 0.0f}, {-1.0f, 0.0f, 0.0f}},
	{{-1.0f,  1.0f, -1.0f}, {0.0f, 1.0f}, {-1.0f, 0.0f, 0.0f}},
	{{-1.0f,  1.0f,  1.0f}, {1.0f, 1.0f}, {-1.0f, 0.0f, 0.0f}},

	// Right face
	{{ 1.0f, -1.0f,  1.0f}, {0.0f, 0.0f}, {1.0f, 0.0f, 0.0f}},
	{{ 1.0f,  1.0f,  1.0f}, {0.0f, 1.0f}, {1.0f, 0.0f, 0.0f}},
	{{ 1.0f,  1.0f, -1.0f}, {1.0f, 1.0f}, {1.0f, 0.0f, 0.0f}},
	{{ 1.0f,  1.0f, -1.0f}, {1.0f, 1.0f}, {1.0f, 0.0f, 0.0f}},
	{{ 1.0f, -1.0f, -1.0f}, {1.0f, 0.0f}, {1.0f, 0.0f, 0.0f}},
	{{ 1.0f, -1.0f,  1.0f}, {0.0f, 0.0f}, {1.0f, 0.0f, 0.0f}},

	// Top face
	{{-1.0f,  1.0f, -1.0f}, {0.0f, 0.0f}, {0.0f, 1.0f, 0.0f}},
	{{ 1.0f,  1.0f, -1.0f}, {1.0f, 0.0f}, {0.0f, 1.0f, 0.0f}},
	{{ 1.0f,  1.0f,  1.0f}, {1.0f, 1.0f}, {0.0f, 1.0f, 0.0f}},
	{{ 1.0f,  1.0f,  1.0f}, {1.0f, 1.0f}, {0.0f, 1.0f, 0.0f}},
	{{-1.0f,  1.0f,  1.0f}, {0.0f, 1.0f}, {0.0f, 1.0f, 0.0f}},
	{{-1.0f,  1.0f, -1.0f}, {0.0f, 0.0f}, {0.0f, 1.0f, 0.0f}},

	// Bottom face
	{{-1.0f, -1.0f, -1.0f}, {0.0f, 1.0f}, {0.0f, -1.0f, 0.0f}},
	{{-1.0f, -1.0f,  1.0f}, {0.0f, 0.0f}, {0.0f, -1.0f, 0.0f}},
	{{ 1.0f, -1.0f,  1.0f}, {1.0f, 0.0f}, {0.0f, -1.0f, 0.0f}},
	{{ 1.0f, -1.0f,  1.0f}, {1.0f, 0.0f}, {0.0f, -1.0f, 0.0f}},
	{{ 1.0f, -1.0f, -1.0f}, {1.0f, 1.0f}, {0.0f, -1.0f, 0.0f}},
	{{-1.0f, -1.0f, -1.0f}, {0.0f, 1.0f}, {0.0f, -1.0f, 0.0f}}
	};

	sBuffer = std::make_shared<buffer>(sVertices, drawFreq::staticDraw);
}

void skyBox::draw(GLuint shaderID) {
	sBuffer->draw(drawID::arrayDraw, drawType::triDraw, glm::vec3(0, 0, 0), glm::vec3(0, 0, 0), glm::vec3(1, 1, 1), shaderID);
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