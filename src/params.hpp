#ifndef PARAMS_H
#define PARAMS_H

struct Params {
    char* fname;   /* File name          */
    int   nframes; /* Number of frames   */
    int   npframe; /* Steps per frame    */
    int n; 
    float mass; 
    float h;       /* Particle size      */
    float dt;      /* Time step          */
    float rho0;    /* Reference density  */
    float k;       /* Bulk modulus       */
    float mu;      /* Viscosity          */
    float g;       /* Gravity strength   */
    double bin_size; 
    int num_bins;
};

void default_params(Params* params);
int setter(int argc, char** argv, Params* params);

#endif 
