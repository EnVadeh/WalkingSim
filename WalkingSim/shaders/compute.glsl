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
    vec3 padding3;
};

layout(std140, binding = 1) uniform AtmosphereUBO {
    atmosphereParams atm;
};

layout(rgb16f, binding =0) uniform image2D outImage;


vec2 uvtoTransittanceParams(float u, float v){
    float H = sqrt(atm.atmosphereRad * atm.atmosphereRad - atm.earthRad * atm.earthRad);
    float rho = v * H;
    float r = sqrt(rho * rho + atm.earthRad + atm.earthRad);
    float d_min = r - atm.earthRad;
    float d_max = rho + H;
    float d = d_min +  u * (d_max - d_min);
    float mu = d == 0.0f ? 1.0f : (H * H - rho * rho - d * d) / (2.0f * r * d);
    mu = clamp(mu, -1.0f, 1.0f);
    return vec2(r, mu);
}

float rayIntersectSphere(float r, float mu, float radius) {
	float b = 2.0f * r * mu;
	float c = r * r - radius * radius;

	float det = b * b - 4.0f * c;

	if (det < 0.0f) return 0.0f;
	det = sqrt(det);

	float t1 = (-b + det) * 0.5f;
	float t2 = (-b - det) * 0.5f;

	return (t1 >= 0.0f) ? t1 : t2;
}

vec3 computeTransmittance(float r, float mu) {
	const int STEPS = 250;
	float dt = rayIntersectSphere(r, mu, atm.atmosphereRad);
	dt /= float(STEPS);

	vec3 transmittance = vec3(1.0f, 1.0f, 1.0f);
	float t = 0.0f;

	for (int i = 0; i < STEPS; i++) {
		float r_sample = sqrt(r * r + t * t + 2.0f * r * mu * t);
		float height = r_sample - atm.earthRad;
	
		vec3 extinction = atm.betaR * exp(-height / atm.Hr) + atm.betaM * exp(-height / atm.Hm);
		transmittance *= exp(-extinction * dt);
		t += dt;

	}
	return transmittance;
}

void main() {
    ivec2 pixelCoords = ivec2(gl_GlobalInvocationID.xy);
    float r, mu;
    vec2 uv = vec2(pixelCoords)/vec2(64, 256);
    vec2 temp = uvtoTransittanceParams(uv.x, uv.y);
    r = temp.x;
    mu = temp.y;

    vec3 transmittance = computeTransmittance(r, mu);

    // Write the data to the output texture
    //imageStore(outputTexture, pixelCoords, vec4(transmittance, 1.0));
    imageStore(outputTexture, pixelCoords, vec4(1, 0, 1, 1.0));
}