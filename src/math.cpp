#include <iostream>
#include <fstream>
#include <vector> 
#include <omp.h> 
#include <cmath> 
#include <cstring> 
#include <string> 
#include "math.hpp"

void compute_density(int start, int end, Params* p1, std::vector<Bin>& bins){
    float h  = p1->h;
    float h2 = h*h;
    float h8 = ( h2*h2 )*( h2*h2 );
    float C  = 4 * p1->mass / M_PI / h8;
    #pragma omp parallel for
    for (int bin = start; bin < end; ++bin) {
        if(bins[bin].particles.size() == 0)
            continue;
        for(int j=0; j< bins[bin].particles.size(); j++){
            bins[bin].particles[j].rho = 4 * p1->mass / M_PI / h2;
            // Loop through the neighboring bins and print their indices
            for (int neighbor_bin : bins[bin].neighbours) {
                for (Particle& neigh : bins[neighbor_bin].particles) {
                    float dx = bins[bin].particles[j].x- neigh.x;
                    float dy = bins[bin].particles[j].y- neigh.y;
                    float r2 = dx*dx + dy*dy;
                    float z  = h2-r2;
                    if (z > 0 && bins[bin].particles[j].id != neigh.id) {
                        float rho_ij = C*z*z*z;
                        bins[bin].particles[j].rho += rho_ij;
                    }
                }
            }
        }
    }
}

void compute_accel( int start, int end, Params* p1, std::vector<Bin>& bins) {
    // Unpack basic parameters
     float h    = p1->h;
     float rho0 = p1->rho0;
     float k    = p1->k;
     float mu   = p1->mu;
     float g    = p1->g;
     float mass = p1->mass;
     float h2   = h*h;

    // Compute density
    compute_density(start, end, p1, bins);

    for (int bin = 0; bin < (p1->num_bins*p1->num_bins); ++bin) {
        if (bins[bin].particles.size() == 0)
            continue;
        for (int j = 0; j < bins[bin].particles.size(); j++) {
            bins[bin].particles[j].ax=0;
            bins[bin].particles[j].ay=-g;
        }
    }

    // Constants for interaction term
    float C0 = mass / M_PI / ( (h2)*(h2) );
    float Cp =  15*k;
    float Cv = -40*mu;
    
    #pragma omp parallel for
    for (int bin = start; bin < end; ++bin) {
        if (bins[bin].particles.size() == 0)
            continue;

        for (int j = 0; j < bins[bin].particles.size(); j++) {
            float rhoi = bins[bin].particles[j].rho;
            
            // Loop through the neighboring bins and print their indices
            for(int neighbor_bin : bins[bin].neighbours) {
                for (Particle& neigh : bins[neighbor_bin].particles) {
                    float dx = bins[bin].particles[j].x - neigh.x;
                    float dy = bins[bin].particles[j].y - neigh.y;
                    float r2 = dx * dx + dy * dy;
                    
                    if (r2 < h2 && bins[bin].particles[j].id != neigh.id) {
                        const float rhoj = neigh.rho;
                        float q = sqrt(r2) / h;
                        float u = 1 - q;
                        float w0 = C0 * u / (rhoi * rhoj);
                        float wp = w0 * Cp * (rhoi + rhoj - 2 * rho0) * u / q;
                        float wv = w0 * Cv;
                        float dvx = bins[bin].particles[j].vx - neigh.vx;
                        float dvy = bins[bin].particles[j].vy - neigh.vy;
                        // Update the particle's acceleration in the array
                        bins[bin].particles[j].ax += (wp * dx + wv * dvx);
                        bins[bin].particles[j].ay += (wp * dy + wv * dvy);
                    }
                }
            }
        }
    }
}

