#version 430 core

layout(local_size_x = 16, local_size_y = 16) in;

const int TRANSMITTANCE_TEXTURE_WIDTH = 1024; 
const int TRANSMITTANCE_TEXTURE_HEIGHT = 1024;
const double PI = 3.1415;
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

layout(rgba32f, binding = 0) uniform image2D transmittanceLUT;

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

dvec2 densityAtHeight(double height){
    double rayleighDensity = exp(float(-height)/atm.Hr);
    double mieDensity = exp(float(-height)/atm.Hm);
    return dvec2(rayleighDensity, mieDensity);
}

double distanceToBoundary(float H, float theta){
    
         
    return 0;
}

double vFromDir(dvec3 dir){
    double cosTheta = dot(dir, vec3(0, 1, 0)); //convert the direction to -1..+1
    cosTheta = -cosTheta; //prepare for normalisation
    return (cosTheta + 1.0f)/(2.0f); //normalise
}

double thetaFromV(double v){
    double NN = v * 2.0f - 1.0f; //un-normalising the v to go from 0..1 to -1..1 where 1 is straight up (cos 0), -1 is straight down(cos 180)
    NN = -NN; //since the uv is the opposite of -1 to +1 (0 is straight up, 1 is straight down) Look at vFromDir for help
    return acos(float(NN));
}

/*
To figure out the intersection length to the circumference/atmosphere, we can do so through triangle cosine law. //Thank you old man Taniwha from discord for helping me out
think X, Y, Z as sides of a triangle, where X is the segment from center of the concentric circles (atmosphere and earth) to the height presently being worked on
Z is the segment from the height to the atmosphere circumference at an angle theta with vertical linear
and Y is the segment from concentric center to the point of intersection between Z segment and the circumference (equal to the radius of the atmosphere)
Using Cosine law, Y^2 = Z^2 - 2XZcos(pi-theta) + X^2
                  Z^2 + Z * -2Xcos(pi-theta) + (X^2 - Y^2)
Now, using algebraic formula, we solve for Z
                  
*/
double intersectionLength(double alt, double theta){
    double a = 1;
    double X_segment = atm.earthRad + alt;
    double b = - 2 * X_segment * cos(float(PI - theta));
    double c = (X_segment * X_segment - atm.atmosphereRad * atm.atmosphereRad);
    double d = b * b - 4 * a * c;
    if(d < 0.0)
        return -1.0;
    double x1 = (-b + sqrt(d))/(2 * a);
    double x2 = (-b - sqrt(d))/(2 * a);
    if(x1 > 0.0) return x1;
    if(x1 > 0.0) return x2;
    return -1.0;
}

/*
To find maximum angle, we find the tangent to the Earth from that altitude. Remember that the tangent to a circle is always perpendicular to the radius of the circle that intersescts that point of contact. 
Then forms a right angled triangle where the hypoteneuse is the altitude + earthRad, the base is the length from the height to the point of contact with the Earth, and the perpendicular is the radius.
Find the base or 'first segment' of the tangent, and then we can use trigonometry to find the Theta of the right angled triangle, after which we subtract PI/2 by theta and get the desired MaxTheta.

Theta = PI - arcos(sqrt((alt + earthRad)^2 - earthRad^2)/(alt + earthRad))
So, we subtract the theta by a small number since we don't want to touch the planet
*/
double maxTheta(double altitude){ //read the comment above for clarification of variable names
    double hypo = altitude + atm.earthRad;
    double numerator = sqrt(hypo * hypo - atm.earthRad * atm.earthRad); 
    double denom = hypo;
    double theta = PI - acos(float(numerator/denom));
    return (theta - 0.01);
}

dvec3 computeTransmittance(dvec2 frag_coord){
    double alt = (atm.atmosphereRad - atm.earthRad) * frag_coord.x + 0.1;
    double theta = thetaFromV(frag_coord.y);
    double mTheta = maxTheta(alt);
    if(theta > mTheta)
        return vec3(1, 0, 0);
    double maxDist = intersectionLength(alt, theta);
    return dvec3(alt, theta, maxDist);

}

void main() {
    ivec2 pixelCoords = ivec2(gl_GlobalInvocationID.xy);
    ivec2 size = ivec2(TRANSMITTANCE_TEXTURE_WIDTH, TRANSMITTANCE_TEXTURE_HEIGHT);
    if(pixelCoords.x > size.x || pixelCoords.y > size.y){
        return;
    }
    dvec2 frag_coord = dvec2(double(pixelCoords.x)/1023, double(pixelCoords.y)/1023);
    dvec3 transmittance = computeTransmittance(frag_coord);
    // Write the data to the output texture
    imageStore(transmittanceLUT, pixelCoords, vec4(transmittance, 1.0));
}