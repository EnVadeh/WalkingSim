#version 430 core

layout(local_size_x = 16, local_size_y = 16) in;

const int TRANSMITTANCE_TEXTURE_WIDTH = 1024; 
const int TRANSMITTANCE_TEXTURE_HEIGHT = 1024;
const double PI = 3.1415;
const vec3 absorption_extinction = { 7.0e-7, 3.54e-7, 1.52e-7 };
const int NUM_SAMPLES = 50;

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
We only know the Angle at which the ray is travelling, and the height from the center from which the ray starts.
We can make a triangle with one side 'a' as the straight line from the center to the ray Origin, whose value is the original height of the ray Origin (earthRad + altitude), 
another side 'b' is the line form the ray Origin to the current sample point in the theta direction with the vertical, which is the distance from ray origin at altitude alt
to the current sample point distance (step size * number of steps), finally the last side is distance from the center to the current sample point, which would be the 
altitude of the point at the given sample point. Using cosine law with (PI - theta) as the angle between 'a' and 'b', we can calculate the value of 'c'
Using Cosine law, C^2 =  A^2 + B^2 - 2ABcos(pi - theta)
*/
double sampleHeight(double alt, double theta, double sampleDistance){
    double a2 = (atm.earthRad + alt) * (atm.earthRad * alt);
    double b2 = sampleDistance * sampleDistance;
    double c2 = a2 + b2 - 2 * sampleDistance * (atm.earthRad + alt) * cos(float(PI - theta));
    return sqrt(c2);
}

/*
To figure out the intersection length to the circumference/atmosphere, we can do so through triangle cosine law. //Thank you old man Taniwha from discord for helping me out
The Earth and the atmosphere of the Earth can be seen as two concentric circles.
Let's make a triangle with X, Y, Z as the sides, where X is the line segment from center of the concentric circles to the altitude presently being worked on
Z is the line segment from the height to the atmosphere circumference at an angle theta with the vertical line
Y is the line segment from the concentric center to the point of intersection between Z and the circumference of the outer circle(equal to the radius of the atmosphere)
Using Cosine law, Y^2 = Z^2 - 2XZcos(pi-theta) + X^2
                  Z^2 + Z * -2Xcos(pi-theta) + (X^2 - Y^2)
Now, using the quadratic algebraic formula, we solve for Z
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
The maximum vertical angle that a view ray can make from an altitude without the view ray being blocked off by the Earth itself will be infinitismally small degree less than the
vertical angle a tangent makes to the Earth form that altitude.
Remember that the tangent to a circle is always perpendicular to the Radius of the circle. 
Using this, we can form a right angled triangle where the Hypoteneuse is the line segment from the center of the Earth to the current altitdue (earthRad + altitude).
The perpendicular of the right angled triangle will be the line segment from the center to the tangent contact point, with the length of earthRad.
The base of the right angled triangle will then be the 'first segment' of the tangent (the segment of the tangent from the starting point to the point of contact).
Let's assume that the angle of the right angled triangle relative to which we have defined the base and perpendicular is alpha. 
We use pythogoraus theorem to find the value of the base of the triangle. 
We can then use basic trigonometry to find the value of the inner angle.
We know, Cos(alpha) = base/Hypoteneus
                    = segment/(altitude+earthRad)
We need to find maximum theta, which is the outer angle, so it will be PI - alpha
Theta = PI - arcos(sqrt((alt + earthRad)^2 - earthRad^2)/(alt + earthRad))
We then need to subtract a small number to it since we don't want the ray to touch the planet
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
    double stepSize = maxDist/NUM_SAMPLES;
    double rayleigh = 0;
    double mie = 0;

    for(int i = 0; i < NUM_SAMPLES; i++){
        double currentDist = stepSize * i;
        double currentHeight = sampleHeight(alt, theta, currentDist);
        dvec2 density = densityAtHeight(currentHeight);
        density = density * stepSize;
        rayleigh = rayleigh + density.x;
        mie = mie + density.y;
    }

    dvec3 rayleighOutScat = 4 * PI * atm.betaR;
    dvec3 mieOutScat = 4 * PI * atm.betaMext;
    //now this is only for the outscattering to random point, we calculate outscattering ot the sun too and then that's the inscattering afaik
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