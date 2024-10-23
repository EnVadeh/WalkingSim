#version 430 core

struct Light{
	vec4 vLightPos; 
	vec4 vLightDir;
	vec4 vLightColor;
	vec4 vEndPointA;
	vec4 vEndPointB;
};

layout (std140, binding = 0) uniform LightBuffer{
	Light lights[100];
};

in vec2 fTex;

uniform sampler2D water;

out vec4 outColor;

void main(){
	//outColor = vec4(1, 1, 1,1);
	//outColor = texture2D(water, fTex);
	outColor = lights[0].vLightColor;
	//outColor = vec4(fTex, 0, 1);
}