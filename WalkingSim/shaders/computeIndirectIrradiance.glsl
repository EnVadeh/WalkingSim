#version 430 core

layout(local_size_x = 16, local_size_y = 16) in;

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
const int SCATTERING_TEXTURE_WIDTH = SCATTERING_TEXTURE_MU_S_SIZE * SCATTERING_TEXTURE_NU_SIZE;
const int SCATTERING_TEXTURE_HEIGHT = SCATTERING_TEXTURE_MU_SIZE;
const int SCATTERING_TEXTURE_DEPTH = SCATTERING_TEXTURE_R_SIZE;
const ivec3 SCATTERING_SIZE = ivec3(SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT, SCATTERING_TEXTURE_DEPTH);
const float EarthRayleighScaleHeight = 8.0f;
const float EarthMieScaleHeight = 1.2f;
const float miePhaseFunction_g = 0.80;
const float PI = 3.1415;
const float mu_s_min = -0.2076; // I am confusion

uniform int scatteringORDER;

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

layout(rgba32f, binding = 0) uniform image2D transmittanceLUT;
layout(rgba32f, binding = 1) uniform image3D scatteringLUT;
layout(rgba32f, binding = 2) uniform image3D rayleighLUT;
layout(rgba32f, binding = 3) uniform image3D mieLUT;
layout(rgba32f, binding = 4) uniform image2D deltaIrradianceLUT;
layout(rgba32f, binding = 5) uniform image2D scatteringDensityLUT;
layout(rgba32f, binding = 6) uniform image2D irradianceLUT;

float clampCosine(float mu) {
  return clamp(mu, float(-1.0), float(1.0));
}

float clampRadius(float r) {
  return clamp(r, atm.earthRad, atm.atmosphereRad);
}

float clampDistance(float d){
  return max(d, 0.0);
}

float safeSqrt(float num){
  return sqrt(max(num, 0.0));
}

float getLayerDensity(int layer, float altitude) {
  float density = dp[layer].exp_term * exp(dp[layer].exp_scale * altitude) +
      dp[layer].linear_term * altitude + dp[layer].constant_term;
  return clamp(density, float(0.0), float(1.0));
}

float getProfileDensity(int layer, float altitude) {
  return altitude < dp[0].width ? getLayerDensity(0, altitude) : getLayerDensity(layer, altitude);

}

float getTextureCoordFromUnitRange(float x, int texture_size){
    return 0.5 / float(texture_size) + x * (1.0 - 1.0 / float(texture_size));
}

float getUnitRangeFromTextureCoord(float u, int texture_size){
    return (u - 0.5 / float(texture_size)) / (1.0 - 1.0 / float(texture_size));
}

float distanceToTopAtmosphereBoundary(
    float r, float mu) {
    float discriminant = r * r * (mu * mu - 1.0) + atm.atmosphereRad * atm.atmosphereRad;
    return clampDistance(-r * mu + sqrt(max(discriminant, 0.0)));
}

float distanceToBottomAtmosphereBoundary(
    float r, float mu) {
    float discriminant = r * r * (mu * mu - 1.0) + atm.earthRad * atm.earthRad;
    return clampDistance(-r * mu - sqrt(max(discriminant, 0.0)));
}

float distanceToNearestAtmosphereBoudnary(float r, float mu, bool ray_r_mu_intersects_ground){
    if(ray_r_mu_intersects_ground)
        return distanceToBottomAtmosphereBoundary(r, mu);
    else
        return distanceToTopAtmosphereBoundary(r, mu);
}

float rayleighPhaseFunction(float nu) {
  float k = 3.0 / (16.0 * PI);
  return k * (1.0 + nu * nu);
}

float miePhaseFunction(float g, float nu) {
  float k = 3.0 / (8.0 * PI) * (1.0 - g * g) / (2.0 + g * g);
  return k * (1.0 + nu * nu) / pow(1.0 + g * g - 2.0 * g * nu, 1.5);
}

vec4 getScatteringTextureUVWZfromRmuMuSNu(float r, float mu, float mu_s, float nu, bool ray_r_mu_intersects_ground){
    float H = sqrt(atm.atmosphereRad * atm.atmosphereRad - atm.earthRad * atm.earthRad);
    float rho = safeSqrt(r * r - atm.earthRad * atm.earthRad);
    float u_r = getTextureCoordFromUnitRange(rho/H, SCATTERING_TEXTURE_R_SIZE);
    float r_mu = r * mu; 
    float discriminant = r_mu * r_mu - r * r + atm.earthRad * atm.earthRad;//discriminant of the intersection of ray r,mu with the ground
    float u_mu;
    if(ray_r_mu_intersects_ground){
        // Distance to the ground for the ray (r,mu), and its minimum and maximum
        // values over all mu - obtained for (r,-1) and (r,mu_horizon).
        float d = -r_mu - safeSqrt(discriminant);
        float d_min = r - atm.earthRad;
        float d_max = rho;
        u_mu = 0.5 - 0.5 * getTextureCoordFromUnitRange(d_max == d_min ? 0.0 : (d - d_min) / (d_max - d_min), SCATTERING_TEXTURE_MU_SIZE / 2);
    }
    else{
        // Distance to the top atmosphere boundary for the ray (r,mu), and its
        // minimum and maximum values over all mu - obtained for (r,1) and
        // (r,mu_horizon).
        float d = -r_mu + safeSqrt(discriminant + H * H);
        float d_min = atm.atmosphereRad - r;
        float d_max = rho + H;
        u_mu = 0.5 + 0.5 * getTextureCoordFromUnitRange((d - d_min) / (d_max - d_min), SCATTERING_TEXTURE_MU_SIZE / 2);
    }
    float d = distanceToTopAtmosphereBoundary(atm.earthRad, mu_s);
    float d_min = atm.atmosphereRad - atm.earthRad;
    float d_max = H;
    float a = (d - d_min) / (d_max - d_min);
    float A = -2.0 * mu_s_min * atm.earthRad / (d_max - d_min);
    float u_mu_s = getTextureCoordFromUnitRange(max(1.0 - a / A, 0.0) / (1.0 + a), SCATTERING_TEXTURE_MU_S_SIZE);

    float u_nu = (nu + 1.0) / 2.0;
    return vec4(u_nu, u_mu_s, u_mu, u_r);
}

vec3 getScatteringRayleigh(float r, float mu, float mu_s, float nu, bool ray_r_mu_intersects_ground) {
  vec4 uvwz = getScatteringTextureUVWZfromRmuMuSNu(r, mu, mu_s, nu, ray_r_mu_intersects_ground);
  float tex_coord_x = uvwz.x * float(SCATTERING_TEXTURE_NU_SIZE - 1);
  float tex_x = floor(tex_coord_x);
  float lerp = tex_coord_x - tex_x;
  vec3 uvw0 = vec3((tex_x + uvwz.y) / float(SCATTERING_TEXTURE_NU_SIZE),
      uvwz.z, uvwz.w);
  vec3 uvw1 = vec3((tex_x + 1.0 + uvwz.y) / float(SCATTERING_TEXTURE_NU_SIZE),
      uvwz.z, uvwz.w);
  ivec3 iImageCoord0 = ivec3(uvw0 * ivec3(SCATTERING_SIZE));
  ivec3 iImageCoord1 = ivec3(uvw1 * ivec3(SCATTERING_SIZE));
  return vec3(imageLoad(rayleighLUT, iImageCoord0).rgb * (1.0 - lerp) + imageLoad(rayleighLUT, iImageCoord1).rgb * lerp);
}

vec3 getScatteringDelta(float r, float mu, float mu_s, float nu, bool ray_r_mu_intersects_ground) {
  vec4 uvwz = getScatteringTextureUVWZfromRmuMuSNu(r, mu, mu_s, nu, ray_r_mu_intersects_ground);
  float tex_coord_x = uvwz.x * float(SCATTERING_TEXTURE_NU_SIZE - 1);
  float tex_x = floor(tex_coord_x);
  float lerp = tex_coord_x - tex_x;
  vec3 uvw0 = vec3((tex_x + uvwz.y) / float(SCATTERING_TEXTURE_NU_SIZE),
      uvwz.z, uvwz.w);
  vec3 uvw1 = ivec3((tex_x + 1.0 + uvwz.y) / float(SCATTERING_TEXTURE_NU_SIZE),
      uvwz.z, uvwz.w);
  ivec3 iImageCoord0 = ivec3(uvw0 * ivec3(SCATTERING_SIZE));
  ivec3 iImageCoord1 = ivec3(uvw1 * ivec3(SCATTERING_SIZE));
  return vec3(imageLoad(scatteringLUT, iImageCoord0).rgb * (1.0 - lerp) + imageLoad(scatteringLUT, iImageCoord1).rgb * lerp);
}

vec3 getScatteringMie(float r, float mu, float mu_s, float nu, bool ray_r_mu_intersects_ground) {
  vec4 uvwz = getScatteringTextureUVWZfromRmuMuSNu(r, mu, mu_s, nu, ray_r_mu_intersects_ground);
  float tex_coord_x = uvwz.x * float(SCATTERING_TEXTURE_NU_SIZE - 1);
  float tex_x = floor(tex_coord_x);
  float lerp = tex_coord_x - tex_x;
 vec3 uvw0 = vec3((tex_x + uvwz.y) / float(SCATTERING_TEXTURE_NU_SIZE),
      uvwz.z, uvwz.w);
  vec3 uvw1 = vec3((tex_x + 1.0 + uvwz.y) / float(SCATTERING_TEXTURE_NU_SIZE),
      uvwz.z, uvwz.w);
  ivec3 iImageCoord0 = ivec3(uvw0 * ivec3(SCATTERING_SIZE));
  ivec3 iImageCoord1 = ivec3(uvw1 * ivec3(SCATTERING_SIZE));
  return vec3(imageLoad(mieLUT, iImageCoord0).rgb * (1.0 - lerp) + imageLoad(mieLUT, iImageCoord1).rgb * lerp);
}

vec3 getScattering(
    float r, float mu, float mu_s, float nu,
    bool ray_r_mu_intersects_ground,
    int scattering_order) {

  if (scattering_order == 1) {
    vec3 rayleigh = getScatteringRayleigh(r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    vec3 mie = getScatteringMie(r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    return rayleigh * rayleighPhaseFunction(nu) +
        mie * miePhaseFunction(miePhaseFunction_g, nu);
  } else {
    return getScatteringDelta(r, mu, mu_s, nu,
        ray_r_mu_intersects_ground); //turns out delta_multiple_scattering_texture = delta_rayleigh_texture  
  }
}

vec3 computeIndirectIrradiance(float r, float mu_s, int scattering_order) {
  const int SAMPLE_COUNT = 32; 
  const float dphi = PI / float(SAMPLE_COUNT);
  const float dtheta = PI / float(SAMPLE_COUNT);
  vec3 result = vec3(0.0, 0.0, 0.0);
  vec3 omega_s = vec3(sqrt(1.0 - mu_s * mu_s), 0.0, mu_s);
  for (int iphi = 0; iphi < SAMPLE_COUNT / 2; ++iphi) {
    float phi = (float(iphi) + 0.5) * dphi;
    for (int itheta = 0; itheta < 2 * SAMPLE_COUNT; ++itheta) {
      float theta = (float(itheta) + 0.5) * dtheta;
      vec3 omega = vec3(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
      float domega = (dtheta) * (dphi) * sin(theta);
      float nu = dot(omega, omega_s);
      float dw = dtheta * dphi * sin(theta);
      result += getScattering(r, omega.z, mu_s, nu, false /* ray_r_theta_intersects_ground */,scattering_order) * dw * omega.z;
    }
  }
  return result;
}

void getRMuSFromIrradianceTextureUV(vec2 uv, out float r, out float mu_s) {
  float x_mu_s = getUnitRangeFromTextureCoord(uv.x, IRRADIANCE_TEXTURE_WIDTH);
  float x_r = getUnitRangeFromTextureCoord(uv.y, IRRADIANCE_TEXTURE_HEIGHT);
  r = atm.earthRad + x_r * (atm.atmosphereRad - atm.earthRad);
  mu_s = clampCosine(2.0 * x_mu_s - 1.0);
}

vec3 computeIndirectIrradianceTexture(vec2 uv, int ORDER){
    float r;
    float mu_s;
    getRMuSFromIrradianceTextureUV(uv, r, mu_s);
    return computeIndirectIrradiance(r, mu_s, ORDER);
}

void main() {
    ivec2 pixelCoords = ivec2(gl_GlobalInvocationID.xyz);
    ivec2 size = ivec2(IRRADIANCE_TEXTURE_WIDTH, IRRADIANCE_TEXTURE_HEIGHT);
    vec2 frag_coord = vec2(pixelCoords)/vec2(size);
    vec3 irradiance = computeIndirectIrradianceTexture(frag_coord, scatteringORDER);
    imageStore(deltaIrradianceLUT, pixelCoords, vec4(irradiance, 1));
    imageStore(irradianceLUT, pixelCoords, vec4(irradiance, 1));
    
}