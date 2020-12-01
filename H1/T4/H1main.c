// Includes from C lib
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// Our own headers
#include "H1lattice.h"
#include "H1potential.h"
#include "write_to_file.h"

// Constants
#define NDIM 3
#define KB 8.617e-5
#define TO_GPA 160.2176

int main(){
    // Initializing 
    // Time
    double T    = 3.0;
    double dt   = 3.0e-3;
    int n_timesteps = T/dt;

    int save_every = 10;
    int length_saved = n_timesteps/save_every + 1;
    double time[length_saved];

    for(int t = 0; t < length_saved; t++){
        time[t] = t*dt*save_every;
    }

    // Aluminium constants 
    double m_al     = 26.0/9649.0;       // Weight 26u 
    int N           = 4;                 // Unit cells
    int n_atoms     = 4*N*N*N;              

    double pos[n_atoms][NDIM];
    double v[n_atoms][NDIM]; 
    initialize_atoms(pos, n_atoms, "../T4_eq/data/pos_700C_1fs_1.csv");
    initialize_atoms(v, n_atoms, "../T4_eq/data/v_700C_1fs_1.csv");
//    initialize_atoms(pos, n_atoms, "pos_after.csv");
//    initialize_atoms(v, n_atoms, "v_after.csv");


    double V = 0.0;
    get_volume(&V, "../T4_eq/data/final_v_700C_1fs_1.csv");
    double L = cbrt(V);

    // Energy
    double potential_energy[length_saved];
    double kinetic_energy[length_saved];
    double virial[length_saved];

    potential_energy[0] = get_energy_AL(pos, L, n_atoms);
    kinetic_energy[0] = get_kinetic_energy(v, n_atoms, m_al);
    virial[0] = get_virial_AL(pos, L, n_atoms);

    // Forces for verlet
    double f[n_atoms][NDIM];
    for(int i = 0; i < n_atoms; i++){
        for(int j = 0; j < NDIM; j++){
            f[i][j] = 0.0;
        }
    }

    double temperature[length_saved];
    double pressure[length_saved];

    // Number save all
    int all_dim = NDIM*n_atoms;
    double q_trail[length_saved][all_dim];
    double v_trail[length_saved][all_dim];

    for (int i = 0; i < n_atoms; i++){
        q_trail[0][3*i + 0] = pos[i][0];
        q_trail[0][3*i + 1] = pos[i][1];
        q_trail[0][3*i + 2] = pos[i][2];
        v_trail[0][3*i + 0] = v[i][0];
        v_trail[0][3*i + 1] = v[i][1];
        v_trail[0][3*i + 2] = v[i][2];
    }

    int count = 1;

    // Verlet evolving the system
    get_forces_AL(f, pos, L, n_atoms);
    for(int t = 1; t < n_timesteps + 1; t++){
        verlet_timestep(pos, v, f, n_atoms, dt, m_al, L);

        if (t%save_every == 0){
            kinetic_energy[count] = get_kinetic_energy(v, n_atoms, m_al);
            potential_energy[count] = get_energy_AL(pos, L, n_atoms);
            virial[count] = get_virial_AL(pos, L, n_atoms);

            temperature[count] = (2.0/(3.0*n_atoms))*kinetic_energy[count]/KB; 
            pressure[count] =  (TO_GPA/V)*(n_atoms*KB*temperature[count] + virial[count]);

                
            for(int i = 0; i < n_atoms; i++){
                q_trail[count][3*i + 0] = pos[i][0];
                q_trail[count][3*i + 1] = pos[i][1];
                q_trail[count][3*i + 2] = pos[i][2];
                v_trail[count][3*i + 0] = v[i][0];
                v_trail[count][3*i + 1] = v[i][1];
                v_trail[count][3*i + 2] = v[i][2];
            }
            
            count += 1;
            printf("Saving timestep: %d /%d\n", count, length_saved);
        }
    }



    write_TPV(temperature, pressure, time, length_saved, "TPV");
    write_energy(kinetic_energy, potential_energy, time, length_saved);
    print_pos(pos, n_atoms, "pos_after");
    print_pos(v, n_atoms, "v_after");

    FILE *tq = fopen("q_trail.csv", "w");
    for(int i = 0; i < length_saved; i++){
        for(int j = 0; j < all_dim; j++){
            fprintf(tq, "%f,", q_trail[i][j]);
        }
        fprintf(tq, "%f", time[i]);
        fprintf(tq, "\n");
    }
    fclose(tq);
    FILE *tv = fopen("v_trail.csv", "w");
    for(int i = 0; i < length_saved; i++){
        for(int j = 0; j < all_dim; j++){
            fprintf(tv, "%f,", v_trail[i][j]);
        }
        fprintf(tq, "%f", time[i]);
        fprintf(tv, "\n");
    }
    fclose(tv);
}

