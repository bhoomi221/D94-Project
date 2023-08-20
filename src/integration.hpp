#ifndef INT_H
#define INT_H
#include <vector>
#include "particle.hpp"
#define SCALE 2
#define X_MAX (SCALE)
#define Y_MIN 0
#define Y_MAX (SCALE)
#define X_MIN 0
#define DAMP 0.75


void start(std::vector<Bin>& bins, int start_bin, int end_bin, float dt);
void leapfrog_step(std::vector<Bin>& bins, int start_bin, int end_bin, float dt);
void reflect_bc(std::vector<Bin>& bins, int end);
void damp_reflect(int which, float barrier, Particle& pi);

#endif