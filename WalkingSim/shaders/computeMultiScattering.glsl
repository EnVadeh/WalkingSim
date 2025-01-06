#version 430 core

layout(local_size_x = 16, local_size_y = 16, local_size_z = 4) in;

uniform int scatteringORDER;

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
    densityProfileLayer dp[5];//first 2 are atmosphere layers, then rayleigh, and mie, then empty
};

layout(rgba32f, binding = 0) uniform image2D transmittanceLUT;
layout(rgba32f, binding = 1) uniform image3D scatteringLUT;
layout(rgba32f, binding = 2) uniform image3D rayleighLUT;
layout(rgba32f, binding = 3) uniform image3D mieLUT;
layout(rgba32f, binding = 4) uniform image2D deltaIrradianceLUT;
layout(rgba32f, binding = 5) uniform image3D scatteringDensityLUT;
layout(rgba32f, binding = 6) uniform image2D irradianceLUT;
layout(rgba32f, binding = 7) uniform image3D multiScatteringLUT;

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

bool rayIntersectsGround(float r, float mu) {
  return mu < 0.0 && r * r * (mu * mu - 1.0) + atm.earthRad * atm.earthRad >= 0.0;
}

float rayleighPhaseFunction(float nu) {
  float k = 3.0 / (16.0 * PI);
  return k * (1.0 + nu * nu);
}

float miePhaseFunction(float g, float nu) {
  float k = 3.0 / (8.0 * PI) * (1.0 - g * g) / (2.0 + g * g);
  return k * (1.0 + nu * nu) / pow(1.0 + g * g - 2.0 * g * nu, 1.5);
}

vec2 getTransmittanceTextureUVfromRMu(float r, float mu){ //brunetone's implementation has a figure
    float H = sqrt(atm.atmosphereRad * atm.atmosphereRad - atm.earthRad * atm.earthRad);
    float rho = safeSqrt(r * r - atm.earthRad * atm.earthRad);

    float d = distanceToTopAtmosphereBoundary(r, mu);
    float d_min = atm.atmosphereRad - r;
    float d_max = rho + H;
    float x_mu = (d - d_min) / (d_max - d_min);
    float x_r = rho / H;
    ivec2 size = ivec2(TRANSMITTANCE_TEXTURE_WIDTH, TRANSMITTANCE_TEXTURE_HEIGHT);
    return vec2(getTextureCoordFromUnitRange(x_mu, TRANSMITTANCE_TEXTURE_WIDTH), getTextureCoordFromUnitRange(x_r, TRANSMITTANCE_TEXTURE_HEIGHT));
}

vec3 getTransmittanceToTopAtmosphereBoundary(float r, float mu) {
  vec2 uv = getTransmittanceTextureUVfromRMu(r, mu);
  ivec2 texelCoords = ivec2(uv * ivec2(TRANSMITTANCE_TEXTURE_WIDTH, TRANSMITTANCE_TEXTURE_HEIGHT));
  return vec3(imageLoad(transmittanceLUT, texelCoords).rgb);
}

vec3 getTransmittance(float r, float  mu, float d, bool ray_r_mu_intersects_ground) {
  float r_d = clampRadius(safeSqrt(d * d + 2.0 * r * mu * d + r * r));
  float mu_d = clampCosine((r * mu + d) / r_d);

  if (ray_r_mu_intersects_ground) {
    return min(
        getTransmittanceToTopAtmosphereBoundary(r_d, -mu_d) /
        getTransmittanceToTopAtmosphereBoundary(r, -mu),
        vec3(1.0, 1.0, 1.0));
  } else {
    return min(
        getTransmittanceToTopAtmosphereBoundary(r, mu) /
        getTransmittanceToTopAtmosphereBoundary(r_d, mu_d),
        vec3(1.0, 1.0, 1.0));
  }
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

vec3 getScatteringDensity(float r, float mu, float mu_s, float nu, bool ray_r_mu_intersects_ground) {
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
  return vec3(imageLoad(scatteringDensityLUT, iImageCoord0).rgb * (1.0 - lerp) + imageLoad(scatteringDensityLUT, iImageCoord1).rgb * lerp);
}


void getRMuMuSNuFromScatteringTextureUVWZ(vec4 uvwz, out float r, out float mu, out float mu_s, out float nu, out bool ray_r_mu_intersects_ground){
    float H = sqrt(atm.atmosphereRad * atm.atmosphereRad - atm.earthRad * atm.earthRad);
    float rho = H * getUnitRangeFromTextureCoord(uvwz.w, SCATTERING_TEXTURE_R_SIZE);
    r = sqrt(rho * rho + atm.earthRad * atm.earthRad);

    if(uvwz.z < 0.5){
        // Distance to the ground for the ray r,mu and it's minium and maximum
        // values over all mu - obtained for (r,-1) and (r, mu_horizon) - from which
        // we can recover mu:
        float d_min = r - atm.earthRad;
        float d_max = rho;
        float d = d_min + (d_max - d_min) * getUnitRangeFromTextureCoord(1.0 - 2.0 * uvwz.z, SCATTERING_TEXTURE_MU_SIZE / 2);
        mu = d == 0.0 ? float(-1.0) : clampCosine(-(rho * rho + d * d) / (2.0 * r * d));
        ray_r_mu_intersects_ground = true;
    }
    else{
        // Distance to the top atmosphere boundary for the ray (r,mu), and its
        // minimum and maximum values over all mu - obtained for (r,1) and
        // (r,mu_horizon) - from which we can recover mu:
        float d_min = atm.atmosphereRad - r;
        float d_max = rho + H;
        float d = d_min + (d_max - d_min) * getUnitRangeFromTextureCoord(2.0 * uvwz.z - 1.0, SCATTERING_TEXTURE_MU_SIZE / 2);
        mu = d == 0.0 ? float(1.0) : clampCosine((H * H - rho * rho - d * d) / (2.0 * r * d));
        ray_r_mu_intersects_ground = false;
        }

    float x_mu_s = getUnitRangeFromTextureCoord(uvwz.y, SCATTERING_TEXTURE_MU_S_SIZE);
    float d_min = atm.atmosphereRad - atm.earthRad; //found a mistake
    float d_max = H;
    float A = -2.0 * mu_s_min * atm.earthRad / (d_max - d_min);
    float a = (A - x_mu_s * A) / (1.0 + x_mu_s + A);
    float d = d_min + min(a, A) * (d_max - d_min);
    mu_s = d == 0.0 ? float(1.0) : clampCosine((H * H - d * d) / (2.0 * atm.earthRad * d));
    nu = clampCosine(uvwz.x * 2.0 - 1.0);
}


void getRMuMuSNuFromScatteringTextureFragCoord(vec3 frag_coord, out float r, out float mu, out float mu_s, out float nu, out bool ray_r_mu_intersects_ground){
    const vec4 SCATTERING_TEXTURE_SIZE = vec4(SCATTERING_TEXTURE_NU_SIZE - 1, SCATTERING_TEXTURE_MU_S_SIZE, SCATTERING_TEXTURE_MU_SIZE, SCATTERING_TEXTURE_R_SIZE);
    float frag_coord_nu = floor(frag_coord.x / float(SCATTERING_TEXTURE_MU_S_SIZE));
    float frag_coord_mu_s = mod(frag_coord.x, float(SCATTERING_TEXTURE_MU_S_SIZE));
    vec4 uvwz = vec4(frag_coord_nu, frag_coord_mu_s, frag_coord.y, frag_coord.z) / SCATTERING_TEXTURE_SIZE;
    getRMuMuSNuFromScatteringTextureUVWZ(uvwz, r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    nu = clamp(nu, mu * mu_s - sqrt((1.0 - mu * mu) * (1.0 - mu_s * mu_s)), mu * mu_s + sqrt((1.0 - mu * mu) * (1.0 - mu_s * mu_s)));
}

vec3 integrand(float r, float mu, float mu_s, float nu, float d, bool ray_r_mu_intersects_ground){
    float r_i = sqrt(r * r + d * d + 2.0 * r * mu * d);
    float mu_i = (r * mu + d) / r_i;
    float mus_i = (nu * d + mu_s * r) / r_i;
    return getScatteringDensity(r_i, mu_i, mus_i, nu, ray_r_mu_intersects_ground) * getTransmittance(r, mu, d, ray_r_mu_intersects_ground);
}

vec3 computeMultiScattering(float r, float mu, float mu_s, float nu, int scattering_order, bool ray_r_mu_intersects_ground){
    const int SAMPLE_COUNT = 50;
    float dx = distanceToNearestAtmosphereBoudnary(r, mu, ray_r_mu_intersects_ground) / float(SAMPLE_COUNT);
    float x_i = 0.0;
    vec3 rayleigh_mie_sum = vec3(0, 0, 0);
    vec3 raymie_i = integrand(r, mu, mu_s, nu, 0.0, ray_r_mu_intersects_ground);
    for(int i = 1; i <= SAMPLE_COUNT; ++i){
        float x_j = float (i) * dx;
        vec3 raymie_j = integrand(r, mu, mu_s, nu, x_j, ray_r_mu_intersects_ground);
        rayleigh_mie_sum += (raymie_i + raymie_j) / 2.0f * dx;
        x_i = x_j;
        raymie_i = raymie_j;
    }
    return rayleigh_mie_sum;
}

vec3 computeMultiScatteringTexture(vec3 frag_coord){
    float r;
    float mu;
    float mu_s;
    float nu;
    bool ray_r_mu_intersects_ground;
    getRMuMuSNuFromScatteringTextureFragCoord(frag_coord, r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    return computeMultiScattering(r, mu, mu_s, nu, scatteringORDER, ray_r_mu_intersects_ground);
}

void main() {
    ivec3 pixelCoords = ivec3(gl_GlobalInvocationID.xyz);
    vec3 uv = vec3(pixelCoords) / vec3(SCATTERING_SIZE);
    vec3 multiScattering = computeMultiScatteringTexture(uv);
    imageStore(multiScatteringLUT, pixelCoords, vec4(multiScattering, 1.0f));
}