#ifndef MATH_H
#define MATH_H
#include <iostream>
#include <fstream>
#include <vector> 
#include <cmath> 
#include "params.hpp"
#include "particle.hpp"



void compute_density(int start, int end, Params* p1, std::vector<Bin>& bins);
void compute_accel(int start, int end, Params* p1, std::vector<Bin>& bins);
#endif