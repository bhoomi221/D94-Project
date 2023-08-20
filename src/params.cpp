#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include "params.hpp"


void default_params(Params* params) {
    params->fname   = "run.out";
    params->nframes = 200;
    params->npframe = 200;
    params->dt      = 1e-4f;
    params->h       = 0.1f;
    params->rho0    = 1000;
    params->k       = 1e3f;
    params->mu      = 0.2f;
    params->g       = 9.8f;
    params->n       = 100;
    params->mass    = 1.0f;
}


static void print_usage()
{
    Params param;
    default_params(&param);
    fprintf(stderr, 
            "nbody\n"
            "\t-o: output file name (%s)\n"
            "\t-F: number of frames (%d)\n"
            "\t-f: steps per frame (%d)\n"
            "\t-t: time step (%e)\n"
            "\t-s: particle size (%e)\n"
            "\t-d: reference density (%g)\n"
            "\t-k: bulk modulus (%g)\n"
            "\t-v: dynamic viscosity (%g)\n"
            "\t-g: gravitational strength (%g)\n",
            param.fname, param.nframes, param.npframe,
            param.dt, param.h, param.rho0,
            param.k, param.mu, param.g);
}


int setter(int argc, char** argv, Params* params){

    default_params(params);

    // Vector to store non flag arguments
    //for debugging
    std::vector<std::string> nonOptionArgs;

    // Parsing command-line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

       switch (arg[0]) {
            case '-':
                switch (arg[1]) {
                    case 'h':
                        print_usage(); 
                        break;
                    case 'F':
                        params->nframes = atoi(argv[++i]);
                        break;
                    case 'f':
                        params->npframe = atoi(argv[++i]);
                        break;
                    case 't':
                        params->dt = atof(argv[++i]);
                        break;
                    case 's':
                        params->h = atof(argv[++i]);
                        break;
                    case 'd':
                        params->rho0 = atof(argv[++i]);
                        break;
                    case 'k':
                        params->k = atof(argv[++i]);
                        break;
                    case 'v':
                        params->mu = atof(argv[++i]);
                        break;
                    case 'g':
                        params->g = atof(argv[++i]);
                        break;
                    case 'o':
                        i++;
                        strcpy(params->fname = (char *) malloc(strlen(argv[i])+1), argv[i]);
                        break;  
                    default:
                        std::cerr << "Invalid option: " << arg << std::endl;
                        return 1;
                }
                break;
            default:
                nonOptionArgs.push_back(arg);
                break;
        }
    }

    return 0;
}

