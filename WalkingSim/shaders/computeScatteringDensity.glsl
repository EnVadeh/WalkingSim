#version 430 core

layout(local_size_x = 16, local_size_y = 16, local_size_z = 4) in;

uniform int scatteringORDER;

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
const float miePhaseFunction_g = 0.85; //assymetry parameter for larger areosols
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
    densityProfileLayer dp[5];//first 2 are atmosphere layers, then rayleigh, and mie, then empty
};

layout(rgba16f, binding = 0) uniform image2D transmittanceLUT;
layout(rgba16f, binding = 1) uniform image3D scatteringLUT;
layout(rgba16f, binding = 2) uniform image3D rayleighLUT;
layout(rgba16f, binding = 3) uniform image3D mieLUT;
layout(rgba16f, binding = 4) uniform image2D deltaIrradianceLUT;
layout(rgba16f, binding = 5) uniform image3D scatteringDensityLUT;
layout(rgba16f, binding = 6) uniform image3D multiScatteringLUT;

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
    ivec2 size = imageSize(transmittanceLUT);
    return vec2(getTextureCoordFromUnitRange(x_mu, TRANSMITTANCE_W), getTextureCoordFromUnitRange(x_r, TRANSMITTANCE_H));
}

vec3 getTransmittanceToTopAtmosphereBoundary(float r, float mu) {
  vec2 uv = getTransmittanceTextureUVfromRMu(r, mu);
  ivec2 texelCoords = ivec2(uv * vec2(imageSize(transmittanceLUT)));
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

void getRMuMuSNuFromScatteringTextureFragCoord(vec3 frag_coord, out float r, out float mu, out float mu_s, out float nu, out bool ray_r_mu_intersects_ground){
    const vec4 SCATTERING_TEXTURE_SIZE = vec4(SCATTERING_TEXTURE_NU_SIZE - 1, SCATTERING_TEXTURE_MU_S_SIZE, SCATTERING_TEXTURE_MU_SIZE, SCATTERING_TEXTURE_R_SIZE);
    float frag_coord_nu = floor(frag_coord.x / float(SCATTERING_TEXTURE_MU_S_SIZE));
    float frag_coord_mu_s = mod(frag_coord.x, float(SCATTERING_TEXTURE_MU_S_SIZE));
    vec4 uvwz = vec4(frag_coord_nu, frag_coord_mu_s, frag_coord.y, frag_coord.z) / SCATTERING_TEXTURE_SIZE;
    getRMuMuSNuFromScatteringTextureUVWZ(uvwz, r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    nu = clamp(nu, mu * mu_s - sqrt((1.0 - mu * mu) * (1.0 - mu_s * mu_s)), mu * mu_s + sqrt((1.0 - mu * mu) * (1.0 - mu_s * mu_s)));
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
  ivec3 texelCoords0 = ivec3(uvw0 * vec3(imageSize(rayleighLUT))); //the scatteringLUT will be different...
  ivec3 texelCoords1 = ivec3(uvw1 * vec3(imageSize(rayleighLUT)));
  return vec3(imageLoad(rayleighLUT, texelCoords0).rgb * (1.0 - lerp) +
      imageLoad(rayleighLUT, texelCoords1).rgb * lerp);
}

vec3 getScatteringDelta(float r, float mu, float mu_s, float nu, bool ray_r_mu_intersects_ground) {
  vec4 uvwz = getScatteringTextureUVWZfromRmuMuSNu(r, mu, mu_s, nu, ray_r_mu_intersects_ground);
  float tex_coord_x = uvwz.x * float(SCATTERING_TEXTURE_NU_SIZE - 1);
  float tex_x = floor(tex_coord_x);
  float lerp = tex_coord_x - tex_x;
  vec3 uvw0 = vec3((tex_x + uvwz.y) / float(SCATTERING_TEXTURE_NU_SIZE),
      uvwz.z, uvwz.w);
  vec3 uvw1 = vec3((tex_x + 1.0 + uvwz.y) / float(SCATTERING_TEXTURE_NU_SIZE),
      uvwz.z, uvwz.w);
  ivec3 texelCoords0 = ivec3(uvw0 * vec3(imageSize(rayleighLUT))); //the scatteringLUT will be different...
  ivec3 texelCoords1 = ivec3(uvw1 * vec3(imageSize(rayleighLUT)));
  return vec3(imageLoad(scatteringLUT, texelCoords0).rgb * (1.0 - lerp) +
      imageLoad(scatteringLUT, texelCoords1).rgb * lerp);
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
  ivec3 texelCoords0 = ivec3(uvw0 * vec3(imageSize(mieLUT))); //the scatteringLUT will be different...
  ivec3 texelCoords1 = ivec3(uvw1 * vec3(imageSize(mieLUT)));
  return vec3(imageLoad(mieLUT, texelCoords0).rgb * (1.0 - lerp) +
      imageLoad(mieLUT, texelCoords1).rgb * lerp);
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

vec2 getIrradianceTextureUvFromRMuS(float r, float mu_s) {
  float x_r = (r - atm.earthRad) /
  float (atm.atmosphereRad - atm.earthRad);
  float x_mu_s = mu_s * 0.5 + 0.5;
  return vec2(getTextureCoordFromUnitRange(x_mu_s, IRRADIANCE_TEXTURE_WIDTH),
              getTextureCoordFromUnitRange(x_r, IRRADIANCE_TEXTURE_HEIGHT));
}

vec3 getIrradiance(float r, float mu_s) {
  vec2 uv = getIrradianceTextureUvFromRMuS(r, mu_s);
  ivec2 texelCoords = ivec2(uv * vec2(imageSize(deltaIrradianceLUT)));
  return vec3(imageLoad(deltaIrradianceLUT, texelCoords).rgb);
}

vec3 computeScatteringDensity(float r, float mu, float mu_s, float nu, int scattering_order) {

  // Compute unit direction vectors for the zenith, the view direction omega and
  // and the sun direction omega_s, such that the cosine of the view-zenith
  // angle is mu, the cosine of the sun-zenith angle is mu_s, and the cosine of
  // the view-sun angle is nu. The goal is to simplify computations below.
  vec3 zenith_direction = vec3(0.0, 0.0, 1.0);
  vec3 omega = vec3(sqrt(1.0 - mu * mu), 0.0, mu);
  float sun_dir_x = omega.x == 0.0 ? 0.0 : (nu - mu * mu_s) / omega.x;
  float sun_dir_y = sqrt(max(1.0 - sun_dir_x * sun_dir_x - mu_s * mu_s, 0.0));
  vec3 omega_s = vec3(sun_dir_x, sun_dir_y, mu_s);

  const int SAMPLE_COUNT = 16;
  const float dphi = PI / float(SAMPLE_COUNT);
  const float dtheta = PI / float(SAMPLE_COUNT);
  vec3 rayleigh_mie = vec3(0.0);// RadianceDensitySpectrum(0.0 * watt_per_cubic_meter_per_sr_per_nm);

  // Nested loops for the integral over all the incident directions omega_i.
  for (int l = 0; l < SAMPLE_COUNT; ++l) {
    float theta = (float(l) + 0.5) * dtheta;
    float cos_theta = cos(theta);
    float sin_theta = sin(theta);
    bool ray_r_theta_intersects_ground =
        rayIntersectsGround(r, cos_theta);

    // The distance and transmittance to the ground only depend on theta, so we
    // can compute them in the outer loop for efficiency.
    float distance_to_ground = 0.0;
    vec3 transmittance_to_ground = vec3(0.0, 0.0, 0.0);
    vec3 ground_albedo = vec3(0.0, 0.0, 0.0);
    if (ray_r_theta_intersects_ground) {
      distance_to_ground =
          distanceToBottomAtmosphereBoundary(r, cos_theta);
      transmittance_to_ground =
          getTransmittance(r, cos_theta, distance_to_ground, true /* ray_intersects_ground */);
      ground_albedo = vec3(0.1);//can change this in the future maybe
    }
    for (int samp = 0; samp < 2 * SAMPLE_COUNT; ++samp) { //had to change sample variable to samp because new GLSL
      float  phi = (float(samp) + 0.5) * dphi;
      vec3 omega_i = vec3(cos(phi) * sin_theta, sin(phi) * sin_theta, cos_theta);
      float domega_i = (dtheta) * (dphi) * sin(theta);

      // The radiance L_i arriving from direction omega_i after n-1 bounces is
      // the sum of a term given by the precomputed scattering texture for the
      // (n-1)-th order:
      float  nu1 = dot(omega_s, omega_i);
      vec3 incident_radiance = getScattering(r, omega_i.z, mu_s, nu1, ray_r_theta_intersects_ground, scattering_order - 1);
          ivec3 pixelCoords = ivec3(gl_GlobalInvocationID.xyz);
        imageStore(scatteringDensityLUT, pixelCoords, vec4(incident_radiance, 0));

      // and of the contribution from the light paths with n-1 bounces and whose
      // last bounce is on the ground. This contribution is the product of the
      // transmittance to the ground, the ground albedo, the ground BRDF, and
      // the irradiance received on the ground after n-2 bounces.
      vec3 ground_normal =
          normalize(zenith_direction * r + omega_i * distance_to_ground);
      vec3 ground_irradiance = getIrradiance(atm.earthRad,
          dot(ground_normal, omega_s));
      incident_radiance += transmittance_to_ground *
          ground_albedo * (1.0 / (PI)) * ground_irradiance;

      // The radiance finally scattered from direction omega_i towards direction
      // -omega is the product of the incident radiance, the scattering
      // coefficient, and the phase function for directions omega and omega_i
      // (all this summed over all particle types, i.e. Rayleigh and Mie).
      float nu2 = dot(omega, omega_i);
      float rayleigh_density = getProfileDensity(2, r - atm.earthRad);
      float mie_density = getProfileDensity(3, r - atm.earthRad);
      rayleigh_mie += incident_radiance * (
          atm.betaR * rayleigh_density *
              rayleighPhaseFunction(nu2) +
          atm.betaM * mie_density *
              miePhaseFunction(miePhaseFunction_g, nu2)) *
          domega_i;
    }
  }
  return rayleigh_mie;
}


vec3 computeScatteringDensityTexture(vec3 frag_coord, int scattering_order) {
  float r;
  float mu;
  float mu_s;
  float nu;
  bool ray_r_mu_intersects_ground;
  getRMuMuSNuFromScatteringTextureFragCoord(frag_coord, r, mu, mu_s, nu, ray_r_mu_intersects_ground);
  return computeScatteringDensity(r, mu, mu_s, nu, scattering_order);
}

void main() {
    ivec3 pixelCoords = ivec3(gl_GlobalInvocationID.xyz);
    vec3 frag_coord = vec3(pixelCoords);
    vec3 scattering_density = computeScatteringDensityTexture(vec3(frag_coord.xy, frag_coord.z + 0.5), scatteringORDER); //3 = scattering order, I have to go from 2-4... find an intuitive way to run this computeShader
    //imageStore(scatteringDensityLUT, pixelCoords, vec4(scattering_density, 0));
}