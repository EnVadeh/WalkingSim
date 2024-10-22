#version 430 core

layout(location = 0) in vec3 pos;

uniform mat4 matProjView;
uniform mat4 matModel;

void main(){
	gl_Position = matProjView * matModel * vec4(pos, 1.0);
}