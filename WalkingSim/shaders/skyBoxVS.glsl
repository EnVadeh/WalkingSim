#version 430 core

layout(location = 0) in vec3 pos;
layout(location = 1) in vec2 tex;
layout(location = 2) in vec3 norm;

uniform mat4 matProj;
uniform mat4 matView;

out vec3 fNorm;

void main(){
	fNorm = pos;
	mat4 matProjView = matProj * mat4(mat3(matView));
	gl_Position = (matProjView * vec4(pos, 1.0)).xyww;
	//gl_Position = vec4(pos, 1);
}