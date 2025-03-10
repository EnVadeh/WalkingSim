#version 430 core

layout(location = 0) in vec3 pos;
layout(location = 1) in vec2 tex;
layout(location = 2) in vec3 norm;

out vec2 fTex;

void main() {
	gl_Position = vec4(pos, 1.0);
	fTex = tex;
}