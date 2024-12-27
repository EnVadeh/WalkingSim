#version 430 core

in vec3 fNorm;

uniform sampler2D transmittanceLUT;
uniform sampler3D scatteringLUT;
uniform sampler2D deltaIrradianceLUT;
uniform sampler3D rayleighLUT;
uniform sampler3D mieLUT;
uniform sampler3D scatteringDensityLUT;
uniform vec3 vCamPos;

const float EARTH_RADIUS = 6371e3;
const float ATMOSPHERE_RADIUS = 6471e3;

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

layout(location = 0) out vec3 outColor;


void main(){
	vec2 fragc = gl_FragCoord.xy;
	vec2 uv = fragc/1000.0f;
	outColor = texture(scatteringDensityLUT, vec3(uv, 0.2)).xyz;
}