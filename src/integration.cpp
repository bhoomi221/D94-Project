#include "integration.hpp"
#include <iostream>

void leapfrog_step(std::vector<Bin>& bins, int start_bin, int end_bin, float dt) {
    for (int bin_idx = start_bin; bin_idx < end_bin; ++bin_idx) {
        Bin& bin = bins[bin_idx];
        for (Particle& particle : bin.particles) {
            particle.vhx += particle.ax * dt;
            particle.vhy += particle.ay * dt;
        }
        for (Particle& particle : bin.particles) {
            particle.vx = particle.vhx + particle.ax * dt / 2;
            particle.vy = particle.vhy + particle.ay * dt / 2;
        }
        for (Particle& particle : bin.particles) {
            particle.x += particle.vhx * dt;
            particle.y += particle.vhy * dt;

        }
        
        // Apply boundary conditions to the particles within the bin
    }
    reflect_bc(bins, end_bin);
}



void start(std::vector<Bin>& bins, int start_bin, int end_bin, float dt){
    for (int bin_idx = start_bin; bin_idx < end_bin; ++bin_idx) {
        Bin& bin = bins[bin_idx];
        for (Particle& particle :  bin.particles) {
            particle.vhx = particle.vx + particle.ax * dt / 2;
            particle.vhy = particle.vy + particle.ay * dt / 2;
        }
        for (Particle& particle :  bin.particles) {
            particle.vx += particle.ax  * dt;
            particle.vy += particle.ay  * dt;
        }
        for (Particle& particle :  bin.particles) {
            particle.x += particle.vhx * dt;
            particle.y += particle.vhy * dt;
        }
    }
    reflect_bc(bins, end_bin);
}



void damp_reflect(int which, float barrier, Particle& pi){

    if(which == 0 ){ // vertical bounday case (x)
        // Ignore degenerate cases
        if (pi.vx == 0)
            return;

        // Scale back the distance traveled based on time from collision
        float tbounce = (pi.x-barrier)/pi.vx;
        pi.x -= pi.vx*(1-DAMP)*tbounce;
        pi.y -= pi.vy*(1-DAMP)*tbounce;

        // Reflect the position and velocity
        pi.x  = 2*barrier-pi.x;
        pi.vx  = -pi.vx;
        pi.vhx = -pi.vhx;

        // Damp the velocities
        pi.vx *= DAMP;  pi.vhx *= DAMP;
        pi.vy *= DAMP;  pi.vhy *= DAMP;

    } else { 
        // Ignore degenerate cases
        if (pi.vy == 0)
            return;

        // Scale back the distance traveled based on time from collision
        float tbounce = (pi.y-barrier)/pi.vy;
        pi.x -= pi.vx*(1-DAMP)*tbounce; 
        pi.y -= pi.vy*(1-DAMP)*tbounce;

        // Reflect the position and velocity
        pi.y = 2*barrier-pi.y;
        pi.vy = -pi.vy;
        pi.vhy = -pi.vhy;

        // Damp the velocities
        pi.vx *= DAMP;  pi.vhx *= DAMP;
        pi.vy *= DAMP;  pi.vhy *= DAMP;
    }
}

void reflect_bc(std::vector<Bin>& bins, int end) {
    for (int i = 0; i < end; ++i) {
        Bin& bin = bins[i];
        for(Particle& particle: bin.particles){
            if (particle.x < X_MIN) damp_reflect(0, X_MIN, particle);
            if (particle.x > X_MAX) damp_reflect(0, X_MAX, particle);
            if (particle.y < Y_MIN) damp_reflect(1, Y_MIN, particle);
            if (particle.y > Y_MAX) damp_reflect(1, Y_MAX, particle);
        }
    }
}
