#version 430 core

layout(location = 0) in vec3 pos;
layout(location = 1) in vec2 tex;
layout(location = 2) in vec3 norm;

uniform mat4 matProjView;
uniform mat4 matModel;

layout(rgba32f, binding = 0) uniform image2D transmittanceLUT;

out vec2 fTex;
out vec3 fNorm;
out vec3 vPos;

void main(){
	vec3 tempPos = vec3(pos.x, 0, pos.z);
	gl_Position = matProjView * matModel*  vec4(tempPos, 1.0);
	fTex = tex;
	fNorm = norm;
	vPos = vec3(matModel * vec4(pos, 1.0));
}