#include <iostream>
#include <fstream>
#include <vector> 
#include <cassert> 
#include <cmath>
#include <random>
#include <algorithm> 
#include <mpi.h> 
#include <omp.h>
#include <chrono>
#include "params.hpp"
#include "particle.hpp"
#include "math.hpp"
#include "integration.hpp"

void check_state(std::vector<Bin>& bins, int num_bins)
{
   for (int i =0; i< (num_bins*num_bins); i++){
        for(Particle& p: bins[i].particles){
            float xi = p.x;
            float yi =p.y;
            assert( xi >= 0 || xi <= SCALE );
            assert( yi >= 0 || yi <= SCALE );
        }
    }
}

int circ_indicator(float x, float y)
{
    float dx = (x-(float)(SCALE/2));
    float dy = (y-(float)(SCALE/2));
    float r2 = dx*dx + dy*dy;
    return (r2 < (SCALE/2));
}
int box_indicator(float x, float y)
{
    return (x > 0.1) && (x < (SCALE/2)) && (y < (SCALE/2)) && (y > 0.2);
}
 
int count(Params* p1){
    float h  = p1->h;
    float hh = h/1.2;

    // Count mesh points that fall in indicated region.
    int count = 0;
    for (float x = 0; x < SCALE; x += hh)
        for (float y = 0; y < SCALE; y += hh)
            count += circ_indicator(x,y);

    p1->n =count;
    return count;


}

std::vector<Bin>  place_particles(Params* p1) {
    float h = p1->h;
    float hh = h / 1.2;

    int p = 0;
    std::vector<Bin> bins(p1->num_bins*p1->num_bins);

    for (float x = 0; x < SCALE; x += hh) {
        for (float y = 0; y < SCALE; y += hh) {
            if (circ_indicator(x, y)) {
                Particle temp;
                temp.x= x; 
                temp.y= y; 
                temp.id= p;
                temp.vx=0; 
                temp.vy=0;
                temp.ax=0; 
                temp.ay=0;
                int bin_x = temp.x / p1->bin_size;
                int bin_y = temp.y / p1->bin_size;

                int bin= bin_x + (p1->num_bins*bin_y);
                bins[bin].particles.push_back(temp);  
                p++;
            }

        }
    }
    return bins;
}


void normalize_mass(Params* p1, std::vector<Bin>& bins) {
    float rho2s = 0;
    float rhos = 0;
    compute_density(0, p1->num_bins*p1->num_bins, p1, bins);
    for (int i =0; i< (p1->num_bins*p1->num_bins); i++){
        for(Particle& p: bins[i].particles){
            rho2s += (p.rho) * (p.rho);
            rhos += p.rho;
        }
    }
    p1->mass = ((p1->rho0) * rhos / rho2s);
}

void init_bins(int& num_bins, double& bin_size, float domain_size, float h){
    // Calculate the number of bins 
    num_bins= std::floor(domain_size / (0.2f)); 
    // Calculate bin sizes based on the number of bins
    bin_size = domain_size / num_bins;
}


void temp_populate_bins(std::vector<Particle>& particles, std::vector<Bin>& bins, Params* p1) {
    //delete the old particles 
        for (Bin& bin : bins) {
            bin.particles.clear(); 
            //bin.neighbours.clear();
        }

    // Distribute particles into bins based on their positions
    for (Particle& particle : particles) {
        int bin_x = particle.x / p1->bin_size;
        int bin_y = particle.y / p1->bin_size;
    
        int bin= bin_x + (p1->num_bins*bin_y);
        bins[bin].particles.push_back(particle);   
    }
}

void find_neighbours(int num_bins, std::vector<Bin>& bins) {

    #pragma omp parallel for
    for (int bin_index = 0; bin_index < num_bins* num_bins; ++bin_index) {
        int bin_x = bin_index % num_bins;
        int bin_y = bin_index / num_bins;
        // find neighboring bins for current bin 
        #pragma omp parallel for
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                int neighbor_x = bin_x + dx;
                int neighbor_y = bin_y + dy;
                // check if the neighbor bin is valid 
                if (neighbor_x >= 0 && neighbor_x < num_bins &&
                    neighbor_y >= 0 && neighbor_y < num_bins) {
                    bins[bin_index].neighbours.push_back(neighbor_x + neighbor_y * num_bins);
                }
            }
        }
    }
}


MPI_Datatype createParticleMPIType() {
    MPI_Datatype particle_type;
    int blocklengths[11] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, };
    MPI_Aint displacements_[10] = {offsetof(Particle, id),
                                offsetof(Particle, x),
                                offsetof(Particle, y),
                                offsetof(Particle, vx),
                                offsetof(Particle, vy),
                                offsetof(Particle, vhx),
                                offsetof(Particle, vhy),
                                offsetof(Particle, ax),
                                offsetof(Particle, ay),
                                offsetof(Particle, rho)};
    MPI_Datatype types[10] = {MPI_INT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};
    MPI_Type_create_struct(10, blocklengths, displacements_, types, &particle_type);
    MPI_Type_commit(&particle_type);

    return particle_type;
}

std::vector<Particle> flattenBins(const std::vector<Bin>& bins) {
    std::vector<Particle> flattenedParticles;
    for (const Bin& bin : bins) {
        flattenedParticles.insert(flattenedParticles.end(), bin.particles.begin(), bin.particles.end());
    }
    return flattenedParticles;
}


std::vector<Particle> flattenBinsSp(const std::vector<Bin>& bins, int start, int end) {
    std::vector<Particle> flattenedParticles;
    for (int i =start; i < end; i++ ) {
        flattenedParticles.insert(flattenedParticles.end(), bins[i].particles.begin(), bins[i].particles.end());
    }
    return flattenedParticles;
}
std::vector<Bin> unflattenBins(const std::vector<Particle>& flattenedParticles, const std::vector<int>& particleCounts) {
    std::vector<Bin> bins;
    size_t currentParticleIndex = 0;

    for (int particleCount : particleCounts) {
        Bin bin;
        for (int j = 0; j < particleCount; ++j) {
            bin.particles.push_back(flattenedParticles[currentParticleIndex]);
            ++currentParticleIndex;
        }
        bins.push_back(bin);
    }

    return bins;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    double start_time, end_time;
    if (rank == 0) {
        start_time = MPI_Wtime();
    }
    Params p1;
    if (setter(argc, argv, &p1) != 0)
        exit(-1);
    // Get size and number of bins 
    count(&p1);
    init_bins(p1.num_bins, p1.bin_size, SCALE, p1.h);
    std::vector<Bin> bins(p1.num_bins*p1.num_bins);
    std::vector<Particle> all(p1.n);
    //create datatypes to bcast properly
    MPI_Datatype mpi_particle_type = createParticleMPIType();
    std::vector<int> particle_counts(p1.num_bins*p1.num_bins);
    float mass; 
    if(rank ==0) {
        // Place particles in appropirate bins      
        bins= place_particles(&p1);

        // For each bin create the list of its neighbours 
        find_neighbours(p1.num_bins, bins);
        for(int i=0; i< p1.num_bins*p1.num_bins; i++ ){
            particle_counts[i]=bins[i].particles.size();
        }

        //write the header to the file
        std::ofstream outputFile(p1.fname); 
        outputFile << "SPHView00 " << p1.n << " " << SCALE << std::endl;
        outputFile.close();

        //intialize the mass and compute intial accel for the particles
        normalize_mass(&p1, bins);
        compute_accel(0, p1.num_bins*p1.num_bins, &p1, bins);
        start(bins, 0, p1.num_bins*p1.num_bins, p1.dt);
        //sanity check to make sure all initialized correctly
        check_state(bins, p1.num_bins);
        all = flattenBins(bins); 
        mass= p1.mass;
    } 
    // broadcast the change in mass
    MPI_Bcast(&mass, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    p1.mass= mass;

    // Send the populated bins to other processesby flattening b/c hard to send the struct
    MPI_Bcast(all.data(), all.size(), mpi_particle_type, 0, MPI_COMM_WORLD);
    MPI_Bcast(particle_counts.data(), particle_counts.size(), MPI_INT, 0, MPI_COMM_WORLD);
    bins = unflattenBins(all, particle_counts);
    // Add the neighbors to the Bin (faster look up)
    find_neighbours(p1.num_bins, bins);

    // Determining the valid indices for each process 
    int num_bins_per_process = (p1.num_bins * p1.num_bins) / size;
    std::vector<int> counts(size);
    std::vector<int> displacements(size);
    for (int i = 0; i < size; ++i) {
        counts[i] = num_bins_per_process + (i < (p1.num_bins * p1.num_bins) % size);
        displacements[i] = (i > 0) ? (displacements[i-1] + counts[i-1]) : 0;
    }
    int start_bin = displacements[rank];
    int end_bin = displacements[rank] + counts[rank];

    std::vector<Particle> all_particles(p1.n);
    // Number of frames that will be written to the file
    for(int k =0; k < p1.nframes ; k++){
        //do all math npframe times before writing to the file
        for(int i =0; i< p1.npframe; i++){
            compute_accel(start_bin, end_bin, &p1, bins);
            leapfrog_step(bins, start_bin, end_bin, p1.dt);
            // Flatten to send 
            std::vector<Particle> temp = flattenBinsSp(bins, start_bin, end_bin);
            std::vector<int> recv_counts(size);
            int tsize=temp.size();
            MPI_Allgather(&tsize, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
            std::vector<int> displacements(size);
            int displacement = 0;
            for (int i = 0; i < size; ++i) {
                displacements[i] = displacement;
                displacement += recv_counts[i];
            }
            MPI_Allgatherv(temp.data(), temp.size(), mpi_particle_type, all_particles.data(), recv_counts.data(), displacements.data(), mpi_particle_type, MPI_COMM_WORLD);
            temp_populate_bins(all_particles, bins, &p1);
        } 

        if(rank == 0){
            //intialize the vector that will write to the file 
            std::ofstream outputFile(p1.fname, std::ios::app);
            all.clear();
            for (const Bin& bin : bins) {
                all.insert(all.end(), bin.particles.begin(), bin.particles.end());
            }
            std::sort(all.begin(), all.end(),
                [](const Particle& p1, const Particle& p2) {
                    return p1.id < p2.id;
                });
            for (int i=0; i<all.size(); i++) {
                outputFile << all[i].x << " " << all[i].y << " " << all[i].id << std::endl;
            }
            outputFile.close();
        }
    }
    if(rank == 0){
        end_time = MPI_Wtime();
        std::cout << "Simulation completed in " << end_time - start_time << " seconds.\n";
    }
    // Free the MPI datatype
    MPI_Type_free(&mpi_particle_type);
    MPI_Finalize();
    return 0;
}
