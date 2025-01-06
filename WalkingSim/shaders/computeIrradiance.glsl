#version 430 core

layout(local_size_x = 16, local_size_y = 16) in;

layout(rgba32f, binding = 0) uniform image2D transmittanceLUT;
layout(rgba32f, binding = 4) uniform image2D deltaIrradianceLUT;

const float SUN_ANGULAR_RADIUS = 0.004675; // in radians
const vec3 solar_irradiance = vec3(1.5f);
const int TRANSMITTANCE_W = 256; 
const int TRANSMITTANCE_H = 64;
const int SCATTERING_TEXTURE_R_SIZE = 32;
const int SCATTERING_TEXTURE_MU_SIZE = 128;
const int SCATTERING_TEXTURE_MU_S_SIZE = 32;
const int SCATTERING_TEXTURE_NU_SIZE = 8;
const int IRRADIANCE_TEXTURE_WIDTH = 64;
const int IRRADIANCE_TEXTURE_HEIGHT = 16;
const float EarthRayleighScaleHeight = 8.0f;
const float EarthMieScaleHeight = 1.2f;

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

struct densityProfileLayer { 
  float width;
  float exp_term;
  float exp_scale;
  float linear_term;
  float constant_term;
  float padding1;
  float padding2;
  float padding3;
};

layout(std140, binding = 1) uniform AtmosphereUBO {
    atmosphereParams atm;
};

layout(std140, binding = 2) uniform densityProfileUBO {
    densityProfileLayer dp[4];//first 2 are atmosphere layers, then rayleigh, and mie
};

float clampCosine(float mu) {
  return clamp(mu, float(-1.0), float(1.0));
}

float clampDistance(float d){
  return max(d, 0.0);
}

float getTextureCoordFromUnitRange(float x, int texture_size){
    return 0.5 / float(texture_size) + x * (1.0 - 1.0 / float(texture_size));
}

float getUnitRangeFromTextureCoord(float u, int texture_size){
    return (u - 0.5 / float(texture_size)) / (1.0 - 1.0 / float(texture_size));
}

float distanceToTopAtmosphereBoundary( //don't need distance to bottom boundary for transmittance
    float r, float mu) {
    float discriminant = r * r * (mu * mu - 1.0) + atm.atmosphereRad * atm.atmosphereRad;
    return clampDistance(-r * mu + sqrt(max(discriminant, 0.0)));
}

void GetRMuSFromIrradianceTextureUv(vec2 uv, out float r, out float mu_s) {
  float x_mu_s = getUnitRangeFromTextureCoord(uv.x, IRRADIANCE_TEXTURE_WIDTH);
  float x_r = getUnitRangeFromTextureCoord(uv.y, IRRADIANCE_TEXTURE_HEIGHT);
  r = atm.earthRad + x_r * (atm.atmosphereRad - atm.earthRad);
  mu_s = clampCosine(2.0 * x_mu_s - 1.0);
}

vec2 getTransmittanceTextureUVfromRMu(float r, float mu){ //brunetone's implementation has a figure
    float H = sqrt(atm.atmosphereRad * atm.atmosphereRad - atm.earthRad * atm.earthRad);
    float rho = sqrt(r * r - atm.earthRad * atm.earthRad);
    float d = distanceToTopAtmosphereBoundary(r, mu);
    float d_min = atm.atmosphereRad - r;
    float d_max = rho + H;
    float x_mu = (d - d_min) / (d_max - d_min);
    float x_r = rho / H;
    return vec2(getTextureCoordFromUnitRange(x_mu, TRANSMITTANCE_W), getTextureCoordFromUnitRange(x_r, TRANSMITTANCE_H));
}

vec3 getTransmittanceToTopAtmosphereBoundary(float r, float mu) {
    vec2 uv = getTransmittanceTextureUVfromRMu(r, mu);
    ivec2 texelCoords = ivec2(uv * ivec2(TRANSMITTANCE_W, TRANSMITTANCE_H));
    return imageLoad(transmittanceLUT, texelCoords).rgb;
}

vec3 computeDirectIrradiance(float r, float mu_s) {
    float alpha_s = SUN_ANGULAR_RADIUS;
    // Approximate average of the cosine factor mu_s over the visible fraction of
    // the Sun disc.
    float average_cosine_factor = mu_s < -alpha_s ? 0.0 : (mu_s > alpha_s ? mu_s : (mu_s + alpha_s) * (mu_s + alpha_s) / (4.0 * alpha_s));
    return solar_irradiance * getTransmittanceToTopAtmosphereBoundary(r, mu_s) * average_cosine_factor;
}

vec3 ComputeDirectIrradianceTexture(vec2 frag_coord, vec2 size) {
  float r;
  float mu_s;
  GetRMuSFromIrradianceTextureUv(frag_coord / size, r, mu_s);
  return computeDirectIrradiance(r, mu_s);
}

void main(){
    ivec2 pixelCoords = ivec2(gl_GlobalInvocationID.xy);
    vec2 size = vec2(IRRADIANCE_TEXTURE_WIDTH, IRRADIANCE_TEXTURE_HEIGHT);

    vec2 frag_coord = vec2(pixelCoords);

    vec3 irradiance = ComputeDirectIrradianceTexture(frag_coord, size);
    // Write the data to the output texture
    imageStore(deltaIrradianceLUT, pixelCoords, vec4(irradiance, 1.0));
}