#ifndef PARTICLE_H
#define PARTICLE_H

struct Particle{
    int id;
    float x, y;   //positions 
    float vx, vy; //velocities 
    float vhx, vhy; //half velocities 
    float ax, ay; //accelerations 
    float rho; //density and pressure
};

struct Bin {
    std::vector<Particle> particles;
    std::vector<int> neighbours;
};


#endif