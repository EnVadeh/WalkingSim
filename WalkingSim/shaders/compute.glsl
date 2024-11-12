#version 430 core

layout(local_size_x = 16, local_size_y = 16) in;

layout(rgba16f, binding = 0) uniform image2D outputTexture;

void main() {
    ivec2 pixelCoords = ivec2(gl_GlobalInvocationID.xy);
    
    // Generate some test data
    float r = float(pixelCoords.x) / 256.0;
    float g = float(pixelCoords.y) / 256.0;
    float b = (float(pixelCoords.x) + float(pixelCoords.y)) / 512.0;
    
    // Write the data to the output texture
    imageStore(outputTexture, pixelCoords, vec4(1, 1, 0, 1.0));
}