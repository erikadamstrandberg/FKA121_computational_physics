// Includes from C lib
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_const_mksa.h>

// Our own headers
#include "H1lattice.h"
#include "H1potential.h"
#include "write_to_file.h"

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


// Main program
int main(){
    // Initializing 
    // Time
    double T    = 30;
    double dt   = 1e-3;
    int n_timesteps = T/dt;

    int save_every = 1;
    int length_saved = n_timesteps/save_every + 1;
    double time[length_saved];

    for(int t = 0; t < length_saved; t++){
        time[t] = t*dt*save_every;
    }

    // Aluminium constants 
    double m_al     = 26.0/9649.0;       // Weight 26u 
    double a0       = 4.03;              // Equilbrium constant
    int N           = 4;                 // Unit cells
    int n_atoms     = 4*N*N*N;              

    double pos[n_atoms][NDIM];
    double rand_interval = 0.065*a0; // Random displacement interval
    double L = N*a0;                 // Total length of simulated cube
    double V[length_saved]; 
    V[0] = pow(L, 3);            // Volume of simulated cube
    V[1] = pow(L, 3);

    init_fcc(pos, N, a0);            // Creating pos
    random_displacement(pos, n_atoms, rand_interval); 

    // Energy
    double potential_energy[length_saved];
    double kinetic_energy[length_saved];
    double virial[length_saved];

    potential_energy[0] = get_energy_AL(pos, L, n_atoms);
    kinetic_energy[0] = 0;
    virial[0] = get_virial_AL(pos, L, n_atoms);

    // Forces and velocity for verlet
    double f[n_atoms][NDIM];
    double v[n_atoms][NDIM]; 
    for(int i = 0; i < n_atoms; i++){
        for(int j = 0; j < NDIM; j++){
            f[i][j] = 0.0;
            v[i][j] = 0.0;
        }
    }

    double temperature[length_saved];
    double pressure[length_saved];
    temperature[0] = 0;
    pressure[0] = 0;
    
    // Save timetrail for 5 atoms
    double q1[length_saved][NDIM];
    double q2[length_saved][NDIM];
    double q3[length_saved][NDIM];
    double q4[length_saved][NDIM];
    double q5[length_saved][NDIM];

    // Equilibration
    double start_temp = 2;
    double start_pressure = 4;
    double timestep_temp = start_temp/dt;
    double timestep_pressure = start_pressure/dt;
    
    // Values for temp equilibration
    double T_equil = 500.0;
    double tau_t   = 200.0*dt;
    double alpha_t = 1.0;

    // Values for pressure equilibration
    double P_equil      = 1e-4/TO_GPA;
    double bulk_modulus = 62.0/TO_GPA;
    double kappa_p      = 1.0/bulk_modulus;
    double tau_p        = 400.0*dt;
    double alpha_p      = 1.0;
    

    // Print information before Verlet
    printf("number of timesteps: %d\n", n_timesteps);
    printf("initial potential energy: %f\n", potential_energy[0]);
    printf("initial kinetic energy: %f\n", kinetic_energy[0]);

    printf("number of timesteps: %d\n", n_timesteps);
    printf("initial potential energy: %f\n", potential_energy[0]);
    printf("initial kinetic energy: %f\n", kinetic_energy[0]);


    // Verlet evolving the system
    get_forces_AL(f,pos, L, n_atoms);
    for(int t = 1; t < n_timesteps + 1; t++){
        verlet_timestep(pos, v, f, n_atoms, dt, m_al, L);

        kinetic_energy[t] = get_kinetic_energy(v, n_atoms, m_al);
        potential_energy[t] = get_energy_AL(pos, L, n_atoms);
        virial[t] = get_virial_AL(pos, L, n_atoms);

        temperature[t] = (2.0/(3.0*n_atoms))*kinetic_energy[t]/KB; 
        pressure[t] =  (1.0/V[t])*(n_atoms*KB*temperature[t] + virial[t]);
     
        // Rescaling of the velocites.  
        if((t > timestep_temp && t < timestep_pressure) || (t > 2*timestep_pressure && t < 2.5*timestep_pressure)){
            alpha_t = 1.0 + (2.0*dt/tau_t)*((T_equil - temperature[t])/temperature[t]);
            for(int i = 0; i < n_atoms; i++){
                for(int j = 0; j < NDIM; j++){
                    v[i][j] = sqrt(alpha_t)*v[i][j];
                }
            }
            alpha_p = 1.0;
        }
       
        if((t > timestep_pressure && t < 2*timestep_pressure) || t > 2.5*timestep_pressure){
            alpha_p = 1.0 - kappa_p*(dt/tau_p)*(P_equil - pressure[t]);
            for(int i = 0; i < n_atoms; i++){
                for(int j = 0; j < NDIM; j++){
                    pos[i][j] = cbrt(alpha_p)*pos[i][j];
                }
            }
        }

        V[t+1] = alpha_p*V[t];
        L = cbrt(alpha_p)*L;

        for (int i = 0; i < NDIM; i++){
            q1[t][i] = v[0][i];
            q2[t][i] = v[1][i];
            q3[t][i] = v[2][i];
            q4[t][i] = v[3][i];
            q5[t][i] = v[4][i];
        }

        printf("Saving timestep: %d /%d\n", t, n_timesteps);
    }

    write_TPV(temperature, pressure, V, time, length_saved, "TPV");
    print_pos(pos, n_atoms, "pos_after_equil");
    print_pos(v, n_atoms, "v_after_equil");
    
    FILE *ftrail = fopen("pos_timetrails.csv", "w");
    fprintf(ftrail, "x1,y1,z1,");
    fprintf(ftrail, "x2,y2,z2,");
    fprintf(ftrail, "x3,y3,z3,");
    fprintf(ftrail, "x4,y4,z4,");
    fprintf(ftrail, "x5,y5,z5\n");

    for(int i = 0; i < length_saved; i++){
        fprintf(ftrail, "%f,%f,%f,", q1[i][0], q1[i][1], q1[i][2]);
        fprintf(ftrail, "%f,%f,%f,", q2[i][0], q2[i][1], q2[i][2]);
        fprintf(ftrail, "%f,%f,%f,", q3[i][0], q3[i][1], q3[i][2]);
        fprintf(ftrail, "%f,%f,%f,", q4[i][0], q4[i][1], q4[i][2]);
        fprintf(ftrail, "%f,%f,%f\n", q5[i][0], q5[i][1], q5[i][2]);
    }

    fclose(ftrail);

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
