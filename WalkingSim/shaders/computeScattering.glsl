#version 430 core

layout(local_size_x = 16, local_size_y = 16, local_size_z = 4) in;

const float SUN_ANGULAR_RADIUS = 0.004675; // in radians
const vec3 solar_irradiance = vec3(1.0f);
const int SCATTERING_TEXTURE_R_SIZE = 16;
const int SCATTERING_TEXTURE_MU_SIZE = 16;
const int SCATTERING_TEXTURE_MU_S_SIZE = 16;
const int SCATTERING_TEXTURE_NU_SIZE = 4;
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
    vec3 betaMext;
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

layout(rgba16f, binding = 0) uniform image2D transmittanceLUT;
layout(rgba16f, binding = 1) uniform image3D scatteringLUT;
layout(rgba16f, binding = 2) uniform image3D rayleighLUT;
layout(rgba16f, binding = 3) uniform image3D mieLUT;

float clampCosine(float mu) {
  return clamp(mu, float(-1.0), float(1.0));
}

float clampRadius(float r) {
  return clamp(r, atm.earthRad, atm.atmosphereRad);
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
  //return altitude < dp[layer].width ? getLayerDensity(0, altitude) : getLayerDensity(1, altitude);
    return getLayerDensity(layer, altitude);
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
    return -r * mu + sqrt(max(discriminant, 0.0));
}

float distanceToBottomAtmosphereBoundary(
    float r, float mu) {
    float discriminant = r * r * (mu * mu - 1.0) + atm.earthRad * atm.earthRad;
    return -r * mu - sqrt(max(discriminant, 0.0));
}

float distanceToNearestAtmosphereBoudnary(float r, float mu, bool ray_r_mu_intersects_ground){
    if(ray_r_mu_intersects_ground)
        return distanceToBottomAtmosphereBoundary(r, mu);
    else
        return distanceToTopAtmosphereBoundary(r, mu);
}

vec2 getTransmittanceTextureUVfromRMu(float r, float mu){ //brunetone's implementation has a figure
    float H = sqrt(atm.atmosphereRad * atm.atmosphereRad - atm.earthRad * atm.earthRad);
    float rho = safeSqrt(r * r - atm.earthRad * atm.earthRad);

    float d = distanceToTopAtmosphereBoundary(r, mu);
    float d_min = atm.atmosphereRad - r;
    float d_max = rho + H;
    float x_mu = (d - d_min) / (d_max - d_min);
    float x_r = rho / H;
    ivec2 size = imageSize(transmittanceLUT);
    return vec2(getTextureCoordFromUnitRange(x_mu, size.x), getTextureCoordFromUnitRange(x_r, size.y));
}

vec3 GetTransmittanceToTopAtmosphereBoundary(float r, float mu) {
  vec2 uv = getTransmittanceTextureUVfromRMu(r, mu);
  ivec2 texelCoords = ivec2(uv * vec2(imageSize(transmittanceLUT)));
  return vec3(imageLoad(transmittanceLUT, texelCoords).rgb);
}

vec3 GetTransmittance(float r, float  mu, float d, bool ray_r_mu_intersects_ground) {
  float r_d = clampRadius(safeSqrt(d * d + 2.0 * r * mu * d + r * r));
  float mu_d = clampCosine((r * mu + d) / r_d);

  if (ray_r_mu_intersects_ground) {
    return min(
        GetTransmittanceToTopAtmosphereBoundary(r_d, -mu_d) /
        GetTransmittanceToTopAtmosphereBoundary(r, -mu),
        vec3(1.0, 1.0, 1.0));
  } else {
    return min(
        GetTransmittanceToTopAtmosphereBoundary(r, mu) /
        GetTransmittanceToTopAtmosphereBoundary(r_d, mu_d),
        vec3(1.0, 1.0, 1.0));
  }
}

vec3 GetTransmittanceToSun(float r, float mu_s) {
  float sin_theta_h = atm.earthRad/ r;
  float cos_theta_h = -safeSqrt(max(1.0 - sin_theta_h * sin_theta_h, 0.0));
  return GetTransmittanceToTopAtmosphereBoundary(r, mu_s) *
      smoothstep(-sin_theta_h * SUN_ANGULAR_RADIUS, sin_theta_h * SUN_ANGULAR_RADIUS, mu_s - cos_theta_h);
}

void computeSingleScatteringIntegrand(float r, float mu, float mu_s, float nu, float d, bool ray_r_mu_intersects_ground, out vec3 rayleigh, out vec3 mie){
    float r_d = clampRadius(safeSqrt(d * d + 2.0 * r * mu * d + r * r));
    float mu_s_d = clampCosine((r * mu_s + d * nu) / r_d);
    vec3 transmittance = GetTransmittance(r, mu, d, ray_r_mu_intersects_ground) * GetTransmittanceToSun(r_d, mu_s_d);
    rayleigh = transmittance * getProfileDensity(2, r_d - atm.earthRad);
    mie = transmittance * getProfileDensity(3, r_d - atm.earthRad);
}

void computeSingleScattering(float r, float mu, float mu_s, float nu, bool ray_r_mu_intersects_ground, out vec3 rayleigh, out vec3 mie){
    const int SAMPLE_COUNT = 50;
    float dx = distanceToNearestAtmosphereBoudnary(r, mu, ray_r_mu_intersects_ground) / float(SAMPLE_COUNT);

    vec3 rayleigh_sum = vec3(0.0);
    vec3 mie_sum = vec3(0.0);
    for(int i = 0; i <= SAMPLE_COUNT; ++i){
        float d_i = float(i) * dx;
        vec3 rayleigh_i;
        vec3 mie_i;
        computeSingleScatteringIntegrand(r, mu, mu_s, nu, d_i, ray_r_mu_intersects_ground, rayleigh_i, mie_i);
        float weight_i = (i == 0 || i == SAMPLE_COUNT) ? 0.5 : 1.0;
        rayleigh_sum += rayleigh_i * weight_i;
        mie_sum += mie_i * weight_i;
    }
    rayleigh = rayleigh_sum * dx * solar_irradiance * atm.betaR; 
    mie = mie_sum * dx * solar_irradiance * atm.betaM;
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
    float A = -2.0 * -0.5 * atm.earthRad / (d_max - d_min);
    float u_mu_s = getTextureCoordFromUnitRange(max(1.0 - a / A, 0.0) / (1.0 + a), SCATTERING_TEXTURE_MU_S_SIZE);

    float u_nu = (nu + 1.0) / 2.0;
    return vec4(u_nu, u_mu_s, u_mu, u_r);
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
    float d_min = atm.atmosphereRad - atm.atmosphereRad;
    float d_max = H;
    float A = -2.0 * 0.5 * atm.earthRad / (d_max - d_min);
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

void computeSingleScatteringTexture(vec3 frag_coord, out vec3 rayleigh, out vec3 mie){
    float r;
    float mu;
    float mu_s;
    float nu;
    bool ray_r_mu_intersects_ground;
    getRMuMuSNuFromScatteringTextureFragCoord(frag_coord, r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    computeSingleScattering(r, mu, mu_s, nu, ray_r_mu_intersects_ground, rayleigh, mie);
}

void main() {
    ivec3 pixelCoords = ivec3(gl_GlobalInvocationID.xyz);
    vec3 frag_coord = vec3(pixelCoords);
    vec3 rayleigh;
    vec3 mie;
    computeSingleScatteringTexture(frag_coord, rayleigh, mie); //need layer for 3d coordinates
    vec4 scattering = vec4(rayleigh.rgb, mie.r);
    ivec2 fakeUV = pixelCoords.xy;
    fakeUV.y = imageSize(transmittanceLUT).y - fakeUV.y;
    vec4 transmittance = imageLoad(transmittanceLUT, fakeUV);
    vec4 scatterColor = vec4(0.0);
    scatterColor += vec4(atm.betaR,1.0) * transmittance;
    scatterColor += vec4(atm.betaM, 1.0) * transmittance;
    imageStore(scatteringLUT, pixelCoords, scatterColor);
    imageStore(rayleighLUT, pixelCoords, vec4(rayleigh, 0.0));
    imageStore(mieLUT, pixelCoords, vec4(mie, 0.0));

}