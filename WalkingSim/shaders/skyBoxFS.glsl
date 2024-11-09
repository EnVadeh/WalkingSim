#version 430 core

in vec3 fNorm;
in vec4 gl_FragCoord;

uniform sampler2D transmittanceLUT;

layout(location = 0) out vec3 outColor;

void main(){
	vec2 uv = gl_FragCoord.xy;

	outColor = vec3(texture2D(transmittanceLUT, uv));


}