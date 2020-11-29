// Includes from C lib
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_const_mksa.h>

// Our own headers
#include "H1lattice.h"
#include "H1potential.h"
#include "write_to_file.h"

// Constants
#define NDIM 3
#define KB 8.617e-5
#define TO_GPA 160.2176

void random_displacement(double pos[][NDIM], int n_atoms, double interval);

void verlet_timestep(double pos[][NDIM], double v[][NDIM], double f[][NDIM], 
                     int n_atoms, 
                     double dt, 
                     double m, 
                     double L);

double get_kinetic_energy(double v[][NDIM], int n_atoms, double m);

void write_energy(double *kinetic_energy, double *potential_energy, double *time, int length_saved);

void initialize_atoms(double atoms[][NDIM], int n_atoms, char *filename){
    FILE *fp = fopen(filename, "r");
    char *buf_value = malloc(1000*sizeof(char));
    char *end;

    int i = 0;
    while(fgets(buf_value, n_atoms, fp)){
        atoms[i][0] = strtod(strtok(buf_value,","), &end);
        atoms[i][1] = strtod(strtok(NULL, ","), NULL);
        atoms[i][2] = strtod(strtok(NULL, ","), NULL);
        i += 1;
    }
    fclose(fp);
    free(buf_value);
}

void get_volume(double *V, char *filename){
    FILE *fp = fopen(filename, "r");
    char *buf_value = malloc(100*sizeof(char));
    char *end;

    if (fgets(buf_value, 10, fp) != NULL){    
        *V = strtod(strtok(buf_value, ","), &end);
    }
    fclose(fp);
    free(buf_value);
}

int main(){
    // Initializing 
    // Time
    double T    = 0.1;
    double dt   = 3.0e-3;
    int n_timesteps = T/dt;

    int save_every = 1;
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
    initialize_atoms(pos, n_atoms, "../T3_eq/data/pos_500K_1fs_1.csv");
    initialize_atoms(v, n_atoms, "../T3_eq/data/v_500K_1fs_1.csv");
//    initialize_atoms(pos, n_atoms, "pos_after_equil.csv");
//    initialize_atoms(v, n_atoms, "v_after_equil.csv");


    double V = 0.0;
    get_volume(&V, "../T3_eq/data/final_v_500K_1fs_1.csv");
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

    int count = 0;

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

            count += 1;
                
            for(int i = 0; i < n_atoms; i++){
                q_trail[count][0 + 3*i] = pos[i][0];
                q_trail[count][1 + 3*i] = pos[i][1];
                q_trail[count][2 + 3*i] = pos[i][2];
                v_trail[count][0 + 3*i] = v[i][0];
                v_trail[count][1 + 3*i] = v[i][1];
                v_trail[count][2 + 3*i] = v[i][2];
            }
            printf("Saving timestep: %d /%d\n", t, n_timesteps);
        }
    }

    write_TPV(temperature, pressure, time, length_saved, "TPV");
    write_energy(kinetic_energy, potential_energy, time, length_saved);
    print_pos(pos, n_atoms, "pos_after_equil");
    print_pos(v, n_atoms, "v_after_equil");

    FILE *tq = fopen("q_trail.csv", "w");
    for(int i = 0; i < length_saved; i++){
        for(int j = 0; j < all_dim; j++){
            fprintf(tq, "%f,", q_trail[i][j]);
        }
        fprintf(tq, "\n");
    }
    fclose(tq);
    FILE *tv = fopen("v_trail.csv", "w");
    for(int i = 0; i < length_saved; i++){
        for(int j = 0; j < all_dim; j++){
            fprintf(tv, "%f,", v_trail[i][j]);
        }
        fprintf(tv, "\n");
    }
    fclose(tv);
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
    for(int t = 0; t < length_saved; t++){
        fprintf(fp, "%f,%f,%f\n", kinetic_energy[t], potential_energy[t], time[t]);
    }
}
