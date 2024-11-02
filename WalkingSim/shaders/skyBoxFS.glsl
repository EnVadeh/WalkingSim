#version 430 core

in vec3 fNorm;

layout(location = 0) out vec3 outColor;

void main(){
	outColor = fNorm;
}