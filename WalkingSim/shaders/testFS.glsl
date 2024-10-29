#version 430 core

const float PI = 3.14159265359;
const float QUANTIZE = 255; //Number of non black colors

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
in vec3 fNorm;
in vec3 vPos;

uniform sampler2D noise;
uniform vec3 vCamPos;

layout (location = 0) out vec3 outColor;

float remap_tri(float v)
{
    float orig = v * 2.0 - 1.0;
    v = max(-1.0, orig / sqrt(abs(orig)));
    return v - sign(orig) + 0.5;
}

float distribution(float iso, float roughness){
	float a = roughness * roughness;
	float a2 = a*a;
	float iso2 = iso*iso;
	float D = iso2 * (a2 - 1.0) + 1.0;
	D = PI * D * D;
	return a2/max(D, 0.001);
}

float geoSchlick(float NdotV, float roughness){
	float alpha = roughness * roughness;
	float k = (alpha)/2.0;
	float D = NdotV * (1.0 - k) + k;
	return max(NdotV, 0.001)/D;
}

float geoSmith(vec3 N, vec3 V, vec3 L, float roughness){
	float NdotV = max(dot(N, V), 0.0);
	float NdotL = max(dot(N, L), 0.0);
	float ggx1 = geoSchlick(NdotV, roughness);
	float ggx2 = geoSchlick(NdotL, roughness);
	return ggx1 * ggx2;
}

vec3 fres(float cosTheta, vec3 F){
	return F + (1.0 - F) * pow(1.0 - cosTheta, 5.0); //how muych surface actually refelcts
}

void main(){
	float met = 0.7;
	float rough = 0.7;
	vec3 N = fNorm;

	vec3 lightDir = normalize(vec3(lights[0].vLightDir));
	vec3 camDir = normalize(vCamPos-vPos);
	float fAmbient = 0.3;
	vec2 uvNoise = fTex * (vec2(1000) / vec2(textureSize(noise, 0)));
	vec3 dither = vec3(texture2D(noise, uvNoise));
	vec3 color = vec3(1.0);
	vec3 albedoAmbient = color * fAmbient;
	float reflectance = 1;

	vec3 f0 = vec3(0.3 * reflectance * reflectance);
	f0 = mix(vec3(0.03), color, met);

	vec3 halfDir = normalize(lightDir + camDir);
	float attenuation = 1; //For directional light

	vec3 radiance = vec3(2.0) * attenuation; //give random flux per unit solid angle
	float fresnel = dot(halfDir, lightDir);
	float isoMicro = max(dot(halfDir, N), 0.0);
	
	float NDF = distribution(isoMicro, rough);
	float G = geoSmith(N, camDir, lightDir, rough);
	vec3 F = fres(max(dot(halfDir, camDir), 0.0), f0);

	vec3 nominator = NDF * G * F;

	float denominator = 4.0 * max(dot(N, camDir), 0.0) * max(dot(N, lightDir), 0.0); 

	vec3 specular = nominator/max(denominator, 0.001);
	vec3 kS = F;
	vec3 kD = vec3(1.0) - kS;
	kD *= 1.0 - met;

	float NdotL = max(dot(N, lightDir), 0.0);
	vec3 Lo = (kD * color/ PI + specular) * NdotL * radiance;
	Lo += remap_tri(dither.r)/255;
	Lo = clamp(Lo, 0.0, 1.0);
	outColor = Lo;
	//outColor = vec4(fTex, 0, 1);
}