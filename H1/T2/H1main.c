/*
 H1main.c
 
 Created by Anders Lindman on 2013-10-31.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "H1lattice.h"
#include "H1potential.h"
#include "write_to_file.h"

#define NDIM 3

void random_displacement(double pos[][NDIM], int n_atoms, double interval);

void verlet_timestep(double pos[][NDIM], double v[][NDIM], double f[][NDIM], 
                     int n_atoms, 
                     double dt, 
                     double m, 
                     double L);

double get_kinetic_energy(double v[][NDIM], int n_atoms, double m);

void write_energy(double *kinetic_energy, double *potential_energy, double *time, int length_saved);

/* Main program */
int main(){
    // Initializing 
    // Time
    double T = 100;
    double dt = 0.001;
    int n_timesteps = T/dt;

    int save_every = 10;
    int length_saved = n_timesteps/save_every + 1;
    double time[length_saved];

    for(int t = 0; t < length_saved; t++){
        time[t] = t*dt*save_every;
    }

    // Aluminium constants 
    double m_al = 26.0/9649.0;       // Weight 26u 
    double a0 = 4.03;                // Equilbrium constant
    int N = 4;                       // Unit cells
    int n_atoms = 4*N*N*N;              

    double pos[n_atoms][NDIM];
    double rand_interval = 0.065*a0; // Random displacement interval
    double L = N*a0;                 // Total length of simulated cube
    init_fcc(pos, N, a0);            // Creating pos
    random_displacement(pos, n_atoms, rand_interval); 

    print_pos(pos, n_atoms, "initial_random_displacement");

    // Energy
    double potential_energy[length_saved];
    double kinetic_energy[length_saved];
    potential_energy[0] = get_energy_AL(pos, L, n_atoms);
    kinetic_energy[0] = 0;

    // Forces and velocity for verlet
    double f[n_atoms][NDIM];
    double v[n_atoms][NDIM]; 
    for(int i = 0; i < n_atoms; i++){
        for(int j = 0; j < NDIM; j++){
            f[i][j] = 0.0;
            v[i][j] = 0.0;
        }
    }

    // Print information before Verlet
    printf("number of timesteps: %d\n", n_timesteps);
    printf("initial potential energy: %f\n", potential_energy[0]);
    printf("initial kinetic energy: %f\n", kinetic_energy[0]);

    // Verlet evolving the system
    get_forces_AL(f,pos, L, n_atoms);
    int count = 1;
    for(int t = 1; t < n_timesteps; t++){
        verlet_timestep(pos, v, f, n_atoms, dt, m_al, L);
        if(t%save_every == 0){
            printf("Saving timestep: %d\n", t);
            kinetic_energy[count] = get_kinetic_energy(v, n_atoms, m_al);
            potential_energy[count] = get_energy_AL(pos, L, n_atoms);
            count += 1;
        }
    }

    write_energy(kinetic_energy, potential_energy, time, length_saved); 
    print_pos(pos, n_atoms, "after_verlet");
    /*
     Function that calculates the virial in units of [eV]. pos should be a matrix
     containing the positions of all the atoms, L is the length of the supercell 
     and N is the number of atoms.
    */
    /*
     double virial;
     virial = get_virial_AL(pos, L, N);
    */
    
}

void random_displacement(double pos[][NDIM], int n_atoms, double interval){
    srand(time(NULL));
    double random_value;
    double random_disp;

    for(int i = 0; i < n_atoms; i++){
        for(int j = 0; j < NDIM; j++){
            random_value = (double) rand() / (double) RAND_MAX;
            random_disp = random_value*2*interval-interval;
            pos[i][j] += random_disp;
        }
    }
}

void verlet_timestep(double pos[][NDIM], double v[][NDIM], double f[][NDIM], 
                     int n_atoms, 
                     double dt, 
                     double m, 
                     double L){

    for(int i = 0; i < n_atoms; i++){
        for(int j = 0; j < NDIM; j++){
            v[i][j] += dt*f[i][j]/(2.0*m);
            pos[i][j] += dt*v[i][j];
        }
    }

    get_forces_AL(f,pos, L, n_atoms);

    for(int i = 0; i < n_atoms; i++){
        for(int j = 0; j < NDIM; j++){
            v[i][j] += dt*f[i][j]/(2.0*m);
        }
    }
}

double get_kinetic_energy(double v[][NDIM], int n_atoms, double m){
    double energy = 0.0;
    for(int i = 0; i < n_atoms; i++){
        for(int j = 0; j < NDIM; j++){
            energy += m*pow(v[i][j],2)/2.0;
        }
    }  
    return energy;
}


void write_energy(double *kinetic_energy, double *potential_energy, double *time, int length_saved){
    FILE *fp = fopen("energy.csv", "w");
    fprintf(fp, "kinetic_energy,potential_energy,time\n");
    for(int t = 0; t < length_saved-1; t++){
        fprintf(fp, "%f,%f,%f\n", kinetic_energy[t], potential_energy[t], time[t]);
    }
}
