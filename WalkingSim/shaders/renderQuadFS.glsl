#version 430 core

in vec2 fTex;

uniform sampler2D colorRT;

out vec4 outColor;

void main(){
	outColor = texture2D(colorRT, fTex); //I gotta start dithering in here if I'm doing edge detection, befor ethis pass if not
}