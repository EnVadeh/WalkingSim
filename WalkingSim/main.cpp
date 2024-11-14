#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include "include/utils.hpp"
#include <sstream>
#include <random>
#include "include/camera.hpp"
#include "include/shader.hpp"
#include "include/buffer.hpp"
#include "include/texture.hpp"
#include "include/light.hpp"
#include "include/chunkmanager.hpp"


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

camera* mainCam = nullptr;
void updateDeltaTime() {
	float currentFrame = glfwGetTime();
	deltaTime = currentFrame - lastFrame;
	lastFrame = currentFrame;
}
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
	
	shader computetest("shaders/compute.glsl");
	GLuint computeShader = computetest.createShader();

	shader test("shaders/testVS.glsl", "shaders/testFS.glsl");
	GLuint testShader = test.createShader();
	
	shader RQ("shaders/renderQuadVS.glsl", "shaders/renderQuadFS.glsl");
	GLuint renderProgram = RQ.createShader();

	shader SS("shaders/skyBoxVS.glsl", "shaders/skyBoxFS.glsl");
	GLuint skyProgram = SS.createShader();
	
	textureManager testure;
	testure.loadTexture("E:/NEW_DOanload/grtex.jpg", "water");
	
	textureManager noise;
	noise.loadTexture("E:/NEW_DOanload/WhiteNoiseDithering.png", "noise");
	noise.bindTexture(0, 0, 1, testShader);

	camera myCam(glm::vec3{ 0.0, 1.0, 5.0 }, glm::vec3{ 0.0, 0.0, -1.0 }, -90.0, 0.0);
	mainCam = &myCam;
	//terrain tesst(10, 10);
	lightManager testlights;
	testlights.initLight(glm::vec4(0.7071, 0.7071, 0, 0), glm::vec4(1, 0, 0, 0));
	testlights.turnOn(0);
	testlights.setLights();//binding = 0

	frameBuffer firstpass;
	screenQuad screen(1000, 1000);

	chunkManager cM(myCam);
	atmosphereParams atmosphere;
	auto start = std::chrono::high_resolution_clock::now();
	atmosphereLUTs LUTs(atmosphere); //binding = 1
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
	std::cout << "The time taken for atmosphere LUTs to be claculated: " << duration.count() << std::endl;
	skyBox mSky;

	computeOutput LUT;
	LUT.setup(512, 512);
	glUseProgram(computeShader);
	glDispatchCompute(512 / 16, 512 / 16, 1); //Basically for one texture, z = 1, x = how many groups, y = how many groups
	glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

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
		//tesst.draw(testShader, glm::vec3(0, -1, 0), glm::vec3(1, 1, 1));
		//drawing to the stencil buffer
		glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
		glDepthMask(GL_TRUE);
		glStencilFunc(GL_EQUAL, 1, 0xFF);
		glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
		//drawing stuff to the frame buffer
		glDisable(GL_STENCIL_TEST);
		firstpass.bind();
		glDisable(GL_CULL_FACE);

		noise.bindTexture(0, 0, 1, testShader);
		cM.checkPos(testShader);
		glDepthFunc(GL_LEQUAL);
		//LUT.bind(skyProgram);
		LUT.bind(skyProgram, 0);
		mSky.draw(skyProgram);
		//tesst.draw(testShader, glm::vec3(0, -1.0, 0), glm::vec3(1, 1, 1));
		glDisable(GL_DEPTH_TEST);
		firstpass.sample(renderProgram, 0);
		screen.draw(renderProgram);
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