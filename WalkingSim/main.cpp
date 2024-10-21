#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include "include/utils.hpp"
#include <sstream>
#include <random>
#include "include/camera.hpp"
#include "include/shader.hpp"

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

	glfwSwapInterval(1);
	//put these things in a fucniton so i can just call the draw frame buffers and all these things type shit
	glEnable(GL_DEBUG_OUTPUT);
	while (!glfwWindowShouldClose(window)) {
		GLenum err;
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