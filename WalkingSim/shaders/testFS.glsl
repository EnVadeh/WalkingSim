#version 430 core

uniform sampler2D water;
uniform sampler2D heightmap;

out vec4 outColor;


void main(){
	outColor = vec4(1.0, 1.0,1.0, 0);
}