#version 430 core

layout(local_size_x = 16, local_size_y = 16) in;

const int TRANSMITTANCE_TEXTURE_WIDTH = 256; 
const int TRANSMITTANCE_TEXTURE_HEIGHT = 64;
const vec3 absorption_extinction = { 7.0e-7, 3.54e-7, 1.52e-7 };

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

layout(rgba16f, binding = 0) uniform image2D transmittanceLUT;

layout(std140, binding = 1) uniform AtmosphereUBO {
    atmosphereParams atm;
};

layout(std140, binding = 2) uniform densityProfileUBO {
    densityProfileLayer dp[5];//first 2 are atmosphere layers, then rayleigh, and mie, then empty
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

float ClampCosine(float mu) {
  return clamp(mu, float(-1.0), float(1.0));
}

float distanceToTopAtmosphereBoundary( //don't need distance to bottom boundary for transmittance
    float r, float mu) {
    float discriminant = r * r * (mu * mu - 1.0) + atm.atmosphereRad * atm.atmosphereRad;
    return -r * mu + sqrt(max(discriminant, 0.0));
}

float DistanceToBottomAtmosphereBoundary(
    float r, float mu) {
    float discriminant = r * r * (mu * mu - 1.0) + atm.earthRad * atm.earthRad;
    return -r * mu - sqrt(max(discriminant, 0.0));
}

float getLayerDensity(int layer, float altitude) {
  float density = dp[layer].exp_term * exp(dp[layer].exp_scale * altitude) +
      dp[layer].linear_term * altitude + dp[layer].constant_term;
  return clamp(density, float(0.0), float(1.0));
}

float getProfileDensity(int layer, float altitude) {
    return altitude < dp[0].width ? getLayerDensity(0, altitude) : getLayerDensity(layer, altitude);

}

float ComputeOpticalLengthToTopAtmosphereBoundary(float r, float mu, float scaleHeight) {
    int SAMPLE_COUNT = 500;
    float dx = distanceToTopAtmosphereBoundary(r, mu) / float(SAMPLE_COUNT);
    float result = 0.0f;
    float y_j = exp( -(r - atm.earthRad) / scaleHeight);
    for (int i = 0; i <= SAMPLE_COUNT; ++i) {
        float d_i = float(i) * dx;
    // Distance between the current sample point and the planet center.
        float r_i = sqrt(d_i * d_i + 2.0 * r * mu * d_i + r * r);
    // Number density at the current sample point (divided by the number density
    // at the bottom of the atmosphere, yielding a dimensionless number).
        float y_i = exp(-(r_i - atm.earthRad)/scaleHeight);
    // Sample weight (from the trapezoidal rule).
        //float weight_i = i == 0 || i == SAMPLE_COUNT ? 0.5 : 1.0;
        //result += y_i * weight_i * dx;
        result += (y_j + y_i) / 2.0f * dx;
        y_j = y_i;
  }
  return result;
}


vec3 computeTransmittanceToTopAtmosphereBoundary(float r, float mu){ //basically how much light still is retained when it comes from the sun(source) to the point that we're looking at
    vec3 transmittance = atm.betaR * ComputeOpticalLengthToTopAtmosphereBoundary(r, mu, atm.Hr) + atm.betaM * ComputeOpticalLengthToTopAtmosphereBoundary(r, mu, atm.Hm);
    return exp(-transmittance);
}

/*
vec3 computeTransmittanceToTopAtmosphereBoundary(float r, float mu){ //basically how much light still is retained when it comes from the sun(source) to the point that we're looking at
    return exp(-(atm.betaR * ComputeOpticalLengthToTopAtmosphereBoundary(2, r, mu) + atm.betaMext * ComputeOpticalLengthToTopAtmosphereBoundary(3, r, mu)));
}*/

float getTextureCoordFromUnitRange(float x, int texture_size){
    return 0.5 / float(texture_size) + x * (1.0 - 1.0 / float(texture_size));
}

float getUnitRangeFromTextureCoord(float u, int texture_size){
    return (u - 0.5 / float(texture_size)) / (1.0 - 1.0 / float(texture_size));
}

vec2 getTransmittanceTextureUVfromRMu(float r, float mu){ //brunetone's implementation has a figure
    float H = sqrt(atm.atmosphereRad * atm.atmosphereRad - atm.earthRad * atm.earthRad);
    float rho = sqrt(r * r - atm.earthRad * atm.earthRad);

    float d = distanceToTopAtmosphereBoundary(r, mu);
    float d_min = atm.atmosphereRad - r;
    float d_max = rho + H;
    float x_mu = (d - d_min) / (d_max - d_min);
    float x_r = rho / H;
    ivec2 size = imageSize(transmittanceLUT);
    return vec2(getTextureCoordFromUnitRange(x_mu, size.x), getTextureCoordFromUnitRange(x_r, size.y));

}

void getRMuFromTransmittanceTextureUV(vec2 uv, out float r, out float mu){
    float x_mu = getUnitRangeFromTextureCoord(uv.x,TRANSMITTANCE_TEXTURE_WIDTH); 
    float x_r = getUnitRangeFromTextureCoord(uv.y, TRANSMITTANCE_TEXTURE_HEIGHT); 

    float H = sqrt(atm.atmosphereRad * atm.atmosphereRad - atm.earthRad * atm.earthRad);
    float rho = H * x_r;
    r = sqrt(rho * rho + atm.earthRad * atm.earthRad);
    float d_min = atm.atmosphereRad - r;
    float d_max = rho + H;
    float d = d_min + x_mu * (d_max - d_min);
    mu = d == 0.0 ? 1.0 : (H * H - rho * rho - d * d) / (2.0 * r * d);
    mu = ClampCosine(mu);
}

vec3 computeTransmittanceToTopAtmosphereBoundaryTexture(vec2 frag_coord){
    float r;
    float mu;
    const vec2 TRANSMITTANCE_TEXTURE_SIZE = vec2(TRANSMITTANCE_TEXTURE_WIDTH, TRANSMITTANCE_TEXTURE_HEIGHT);
    getRMuFromTransmittanceTextureUV( frag_coord / TRANSMITTANCE_TEXTURE_SIZE, r, mu);
    return computeTransmittanceToTopAtmosphereBoundary(r, mu);
}

void main() {
    ivec2 pixelCoords = ivec2(gl_GlobalInvocationID.xy);
    ivec2 size = imageSize(transmittanceLUT);

    if(pixelCoords.x >= size.x || pixelCoords.y >= size.y){
        return;
    }
    vec2 frag_coord = vec2(pixelCoords);
    vec3 transmittance = computeTransmittanceToTopAtmosphereBoundaryTexture(frag_coord);
    // Write the data to the output texture
    imageStore(transmittanceLUT, pixelCoords, vec4(transmittance, 1.0));
}