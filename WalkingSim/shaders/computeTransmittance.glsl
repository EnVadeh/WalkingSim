#version 430 core

layout(local_size_x = 16, local_size_y = 16) in;

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

layout(std140, binding = 1) uniform AtmosphereUBO {
    atmosphereParams atm;
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

layout(std140, binding = 2) uniform densityProfileUBO {
    densityProfileLayer dp[2];
};

layout(rgba16f, binding = 0) uniform image2D transmittanceLUT;

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

vec3 computeTransmittanceToTopAtmosphereBoundary(float r, float mu){ //basically how much light still is retained when it comes from the sun(source) to the point that we're looking at
    float dx = distanceToTopAtmosphereBoundary(r, mu);
    int SAMPLE_COUNT = 500;
    float dx_step = dx/ float(SAMPLE_COUNT);
    float x = 0.0;
    vec3 result = vec3(0.0);

    for(int i = 0; i < SAMPLE_COUNT; ++i){
        float r_i = sqrt(r * r + (x + 2.0 * r * mu));
        float dr = dx_step;

        result += atm.betaR * exp(-(r_i - atm.earthRad) / atm.Hr) * dr;
        result += atm.betaMext * exp(-(r_i - atm.earthRad) / atm.Hm) * dr;
        x += dx_step;
    }

    return exp(-result);
}

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
    ivec2 size = imageSize(transmittanceLUT);
    float x_mu = getUnitRangeFromTextureCoord(uv.x, size.x); 
    float x_r = getUnitRangeFromTextureCoord(uv.y, size.y); 

    float H = sqrt(atm.atmosphereRad * atm.atmosphereRad - atm.earthRad * atm.earthRad);
    float rho = H * x_r;
    r = sqrt(rho * rho + atm.earthRad * atm.earthRad);
    float d_min = atm.earthRad - r;
    float d_max = rho + H;
    float d = d_min + x_mu * (d_max - d_min);
    mu = d == 0.0 ? 1.0 : (H * H - rho * rho - d * d) / (2.0 * r * d);
    mu = ClampCosine(mu);
}

vec3 computeTransmittanceToTopAtmosphereBoundaryTexture(vec2 frag_coord, vec2 size){
    float r;
    float mu;
    getRMuFromTransmittanceTextureUV( frag_coord / size, r, mu);
    return computeTransmittanceToTopAtmosphereBoundary(r, mu);
}

void main() {
    ivec2 pixelCoords = ivec2(gl_GlobalInvocationID.xy);
    vec2 size = imageSize(transmittanceLUT);

    if(pixelCoords.x >= size.x || pixelCoords.y >= size.y){
        return;
    }
    vec2 frag_coord = vec2(pixelCoords);

    vec3 transmittance = computeTransmittanceToTopAtmosphereBoundaryTexture(frag_coord, size);
    // Write the data to the output texture
    imageStore(transmittanceLUT, pixelCoords, vec4(transmittance, 1.0));
}