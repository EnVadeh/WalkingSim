#version 430 core

in vec3 fNorm;
in vec3 vPos;

uniform sampler2D transmittanceLUT;
uniform sampler3D scatteringLUT;
uniform sampler3D rayleighLUT;
uniform sampler3D mieLUT;
uniform sampler2D deltaIrradianceLUT;
uniform sampler3D scatteringDensityLUT;
uniform sampler2D irradianceLUT;
uniform sampler3D multiScatteringLUT;
uniform vec3 vCamPos;
uniform mat4 matProj;
uniform mat4 matView;
uniform vec3 vCamDir;
const float SUN_ANGULAR_RADIUS = 0.004675; // in radians
const vec3 solar_irradiance = vec3(1.5f);
const int TRANSMITTANCE_TEXTURE_WIDTH = 256; 
const int TRANSMITTANCE_TEXTURE_HEIGHT = 64;
const int SCATTERING_TEXTURE_R_SIZE = 32;
const int SCATTERING_TEXTURE_MU_SIZE = 128;
const int SCATTERING_TEXTURE_MU_S_SIZE = 32;
const int SCATTERING_TEXTURE_NU_SIZE = 8;
const int IRRADIANCE_TEXTURE_WIDTH = 64;
const int IRRADIANCE_TEXTURE_HEIGHT = 16;
const int SCATTERING_TEXTURE_WIDTH = SCATTERING_TEXTURE_MU_S_SIZE * SCATTERING_TEXTURE_NU_SIZE;
const int SCATTERING_TEXTURE_HEIGHT = SCATTERING_TEXTURE_MU_SIZE;
const int SCATTERING_TEXTURE_DEPTH = SCATTERING_TEXTURE_R_SIZE;
const ivec3 SCATTERING_SIZE = ivec3(SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT, SCATTERING_TEXTURE_DEPTH);
const float EarthRayleighScaleHeight = 8.0f;
const float EarthMieScaleHeight = 1.2f;
const float miePhaseFunction_g = 0.80; //assymetry parameter for larger areosols
const float PI = 3.1415;
const float mu_s_min = -0.2076; // I am confusion
float sun_size = tan(SUN_ANGULAR_RADIUS);
float sun_sizey = cos(SUN_ANGULAR_RADIUS);
vec3 kSphereCenter = vec3(0, 0, 1);
vec3 kSphereAlbedo = vec3(0.8);
vec3 kGroundAlbedo = vec3(0, 0, 0.04);

const float EARTH_RADIUS = 6371e3;
const float ATMOSPHERE_RADIUS = 6471e3;

struct Light{
	vec4 vLightPos; 
	vec4 vLightDir;
	vec4 vLightColor;
	vec4 vEndPointA;
	vec4 vEndPointB;
};

struct atmosphereParams{
    float earthRad;
    float atmosphereRad;
    float Hr;
    float Hm;
    vec3 betaR;
    float padding1;
    vec3 betaM;
    float padding2;
    float meiG;
    float padding3;
    float padding4;
    float padding5;
    vec3 betaMext;
    float padding6;
};

layout(std140, binding = 1) uniform AtmosphereUBO {
    atmosphereParams atm;
};

layout (std140, binding = 0) uniform LightBuffer{
	Light lights[100];
};

layout(location = 0) out vec3 outColor;

const int NUM_SAMPLES = 32;
const vec3 wavelengths = vec3(680.0, 550.0, 440.0);  // RGB wavelengths in nm

vec3 sunDirection = normalize(vec3(lights[0].vLightDir));
//vec3 sunDirection = normalize(vec3(0, -1, 0));
vec3 planetCenter = vec3(0, -6300, 0);
float planetRadius = 6300;

float rayIntersect(vec3 rayOrigin, vec3 rayDir, float sphereRadius) {
    vec3 oc = rayOrigin;
    float b = dot(rayDir, oc);
    float c = dot(oc, oc) - sphereRadius * sphereRadius;
    float discriminant = b * b - c;

    if (discriminant < 0.0) return -1.0; // No intersection
    float t1 = -b - sqrt(discriminant);  // Closest intersection
    float t2 = -b + sqrt(discriminant);  // Farther intersection

    if (t1 > 0.0) return t1; // Return the nearest positive root
    if (t2 > 0.0) return t2; // Return the farther positive root if t1 is invalid
    return -1.0;             // No valid intersection
}

vec2 getDensityRatio(vec3 pos) {
    float height = length(pos - planetCenter) - planetRadius;
    float rayleighDensity = exp(-height / atm.Hr);
    float mieDensity = exp(-height / atm.Hm);
    return vec2(rayleighDensity, mieDensity);
}

vec3 calculateScattering(vec3 start, vec3 dir, float maxDist) {
    float stepSize = maxDist / float(NUM_SAMPLES);
    vec3 rayleighScattering = vec3(0.0);
    vec3 mieScattering = vec3(0.0);
    float opticalDepthR = 0.0;
    float opticalDepthM = 0.0;
    vec2 xy;
    vec3 t;
    for(int i = 0; i < 32; i++) {
        vec3 samplePos = start + dir * (stepSize * float(i));
        vec2 density = getDensityRatio(samplePos);
        opticalDepthR += density.x * stepSize;
        opticalDepthM += density.y * stepSize;
        
        float sunRayLength = rayIntersect(samplePos, sunDirection, 3200);
        vec2 sunDensity = getDensityRatio(samplePos + sunDirection * sunRayLength * 0.5);
        float sunOpticalDepthR = sunDensity.x * sunRayLength;
        float sunOpticalDepthM = sunDensity.y * sunRayLength;
        
        vec3 transmittance = exp(-(atm.betaR * (opticalDepthR + sunOpticalDepthR) +
                                 atm.betaMext * (opticalDepthM + sunOpticalDepthM)));
        rayleighScattering += density.x * transmittance * stepSize;
        mieScattering += density.y * transmittance * stepSize;
        xy = density;
        t = samplePos;
    }

    float cosTheta = dot(dir, sunDirection);
    float rayleighPhase = 3.0 / (16.0 * PI) * (1.0 + cosTheta * cosTheta);
    float miePhase = 3.0 / (8.0 * PI) * ((1.0 - miePhaseFunction_g * miePhaseFunction_g) * (1.0 + cosTheta * cosTheta)) /
                     ((2.0 + miePhaseFunction_g * miePhaseFunction_g) * pow(1.0 + miePhaseFunction_g * miePhaseFunction_g - 2.0 * miePhaseFunction_g * cosTheta, 1.5));
    
    //return vec3(xy,1);
    //return t;
    return (rayleighScattering * rayleighPhase * atm.betaR + mieScattering * miePhase * atm.betaMext) * 20.0;
}

vec3 skyray(vec2 uv, float fieldOfView, float aspectRatio)
{
    float d = 0.5 / tan(fieldOfView / 2.0);
    return vec3((uv.x - 0.5) * aspectRatio, uv.y - 0.5, -d);
}

void main() {
	vec2 fragc = gl_FragCoord.xy;
	vec2 uv = fragc/2048.0f;
	vec3 v_sky_ray = mat3(inverse(matView)) * skyray(uv, 1.0471 , 1.0f);
    vec3 rayDir = normalize(v_sky_ray);
    vec3 rayStart = vec3(vCamPos.x, vCamPos.y + 6300, vCamPos.z);  // Camera position
    float atmosphereIntersect = rayIntersect(rayStart, -rayDir, 6600);
    vec3 scattering = calculateScattering(rayStart, -rayDir, atmosphereIntersect);
    vec3 color = 1.0 - exp(-scattering);  // Tone mapping
    
    //outColor = texture(irradianceLUT, vec2(uv)).xyz + texture(rayleighLUT, vec3(uv, 0.0)).xyz + texture(mieLUT, vec3(uv, 0.0)).xyz + texture(scatteringDensityLUT, vec3(uv, 0.0)).xyz + texture(multiScatteringLUT, vec3(uv, 0.0)).xyz + texture(transmittanceLUT, vec2(uv)).xyz + texture(deltaIrradianceLUT, vec2(uv)).xyz;
    outColor = vec3(color);
    //outColor = texture(transmittanceLUT, vec2(uv)).xyz;
 
 }
 