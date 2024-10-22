#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include "include/utils.hpp"
#include <sstream>
#include <random>
#include "include/camera.hpp"
#include "include/shader.hpp"
#include "include/buffer.hpp"
#include "include/texture.hpp"


void GLAPIENTRY MessageCallback(GLenum source,
	GLenum type,
	GLuint id,
	GLenum severity,
	GLsizei length,
	const GLchar* message,
	const void* userParam)
{
	fprintf(stderr, "GL CALLBACK: %s type = 0x%x, severity = 0x%x, message = %s\n",
		(type == GL_DEBUG_TYPE_ERROR ? "** GL ERROR **" : ""),
		type, severity, message);
}

void updateDeltaTime() {
	float currentFrame = glfwGetTime();
	deltaTime = currentFrame - lastFrame;
	lastFrame = currentFrame;
}

camera* mainCam = nullptr;

void keyCallback(GLFWwindow* window, GLint key, GLint scancode, GLint action, GLint mods) {
	if (action == GLFW_PRESS || action == GLFW_REPEAT) {
		mainCam->processKeyboardInput(key);
	}
}

void mouseCallback(GLFWwindow* window, double xpos, double ypos) {
	static bool firstMouse = true;
	static GLfloat lastX = xpos, lastY = ypos;

	if (firstMouse) {
		lastX = xpos;
		lastY = ypos;
		firstMouse = false;
	}

	GLfloat xoffset = xpos - lastX;
	GLfloat yoffset = lastY - ypos;

	lastX = xpos;
	lastY = ypos;
	mainCam->processMouseInput(xoffset, yoffset);
}

int main() {
	if (!glfwInit()) {
		return -1;
	}
	glfwWindowHint(GLFW_DEPTH_BITS, 24);
	GLFWwindow* window = glfwCreateWindow(1000, 1000, "Window", NULL, NULL);
	if (!window) {
		glfwTerminate();
		return -1;
	}

	glfwMakeContextCurrent(window);
	if (glewInit() != GLEW_OK) {
		return -1;
	}
	
	shader test("shaders/testVS.glsl", "shaders/testFS.glsl");
	GLuint testShader = test.createShader();

	GLfloat testvrs[9] =
	{ -0.5, 0, 0,
	   0, 0.5, 0,
	   0.5, 0, 0
	};

	GLuint testVAO, testVBO;
	glCreateVertexArrays(1, &testVAO);
	glCreateBuffers(1, &testVBO);

	glNamedBufferData(testVBO, sizeof(testvrs), testvrs, GL_STATIC_DRAW);

	glVertexArrayAttribFormat(testVAO, 0, 3, GL_FLOAT, GL_FALSE, 0);
	glEnableVertexArrayAttrib(testVAO, 0);
	glVertexArrayVertexBuffer(testVAO, 0, testVBO, 0, 3 * sizeof(float));

	//std::vector<Vertex> apple;
	//buffer tesss(apple, drawFreq::dynamicDraw);
	//textureManager testure;
	//testure.loadTexture("E:/NEW_DOanload/water.jpg", "water");
	//testure.loadTexture("E:/NEW_DOanload/heightmap.png", "heightmap");
	//testure.bindTexture(0, 0, 2, testShader);
	camera myCam(glm::vec3{ 0.0, 0.0, 0.0 }, glm::vec3{ 0.0, 0.0, -1.0 }, -90.0, 0.0);
	mainCam = &myCam;
	glm::mat4 modeltry = createGeometricToWorldMatrix(glm::vec3(20, 20, -20), glm::vec3(0, 0, 0), glm::vec3(20, 20, 20));
	setUniform(testShader, "matModel", modeltry);
	terrain tesst(10, 10);

	glfwSwapInterval(1);
	//put these things in a fucniton so i can just call the draw frame buffers and all these things type shit
	glEnable(GL_DEBUG_OUTPUT);
	while (!glfwWindowShouldClose(window)) {
		GLenum err;
		updateDeltaTime();
		glfwSetKeyCallback(window, keyCallback);
		glfwSetCursorPosCallback(window, mouseCallback);
		glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
		glEnable(GL_DEPTH_TEST);
		glDepthFunc(GL_LESS);
		glDepthMask(GL_TRUE);
		//shadow buffer shit
		glDisable(GL_DEPTH_TEST);
		glEnable(GL_DEPTH_TEST);
		glDepthMask(GL_TRUE);
		glEnable(GL_CULL_FACE);
		glCullFace(GL_BACK);
		glFrontFace(GL_CCW);
		glEnable(GL_STENCIL_TEST);
		glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
		glDepthMask(GL_FALSE);
		glStencilMask(0xFF);
		glStencilFunc(GL_ALWAYS, 1, 0xFF);
		glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE); 
		//drawing to the stencil buffer
		glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
		glDepthMask(GL_TRUE);
		glStencilFunc(GL_EQUAL, 1, 0xFF);
		glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
		//drawing stuff too the frame buffer
		glDisable(GL_STENCIL_TEST);
		glDisable(GL_CULL_FACE);
		glDisable(GL_DEPTH_TEST);
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glClear(GL_COLOR_BUFFER_BIT);
		glUseProgram(testShader);
		glViewport(0, 0, 1000, 1000);
		glBindVertexArray(testVAO);
		glUseProgram(testShader);
		//glDrawArrays(GL_TRIANGLES, 0, 9);
		tesst.draw(testShader, glm::vec3(0, 0, 0), glm::vec3(20, 20, 20));
		//drawing to the final renderquad
		glfwSwapBuffers(window);
		while ((err = glGetError()) != GL_NO_ERROR)
				printf("OpenGL error: %d\n", err);
		glDebugMessageCallback(MessageCallback, 0);
		glfwPollEvents();
	}
	glfwTerminate();
	return 0;
}