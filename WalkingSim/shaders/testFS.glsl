#version 430 core

in vec2 fTex;

uniform sampler2D water;

out vec4 outColor;

void main(){
	//outColor = vec4(1, 1, 1,1);
	outColor = texture2D(water, fTex);
	//outColor = vec4(fTex, 0, 1);
}