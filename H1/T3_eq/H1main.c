// Includes from C lib
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Our own headers
#include "H1lattice.h"
#include "H1potential.h"
#include "write_to_file.h"

// Constants
// Spatial dimensions
#define NDIM 3

// Boltzmann constant in eV/K
#define KB 8.617e-5

// Pressure from eV/angstrom^3 to GPA
#define TO_GPA 160.2176


// Main program
int main(){
    // Initializing 
    // Time
    double T    = 20;
    double dt   = 1e-3;
    int n_timesteps = T/dt;

    // Decide how many timesteps to save.
    // Not in use during equilbration only during production
    int save_every = 1;
    int length_saved = n_timesteps/save_every + 1;
    double time[length_saved];

    // Generates vector for time
    for(int t = 0; t < length_saved; t++){
        time[t] = t*dt*save_every;
    }

    // Aluminium constants 
    double m_al     = 26.0/9649.0;       // Weight 26u 
    double a0       = 4.03;              // Equilbrium constant
    int N           = 4;                 // Unit cells
    int n_atoms     = 4*N*N*N;              

    // Initializes positions of atoms
    double pos[n_atoms][NDIM];
    double rand_interval = 0.065*a0; // Random displacement interval
    double L = N*a0;                 // Total length of simulated cube
    double V[length_saved]; 
    V[0] = pow(L, 3);            // Volume of simulated cube
    V[1] = pow(L, 3);

    init_fcc(pos, N, a0);            // Creating pos
    random_displacement(pos, n_atoms, rand_interval); 

    // Values for energies
    double potential_energy[length_saved];
    double kinetic_energy[length_saved];
    double virial[length_saved];

    // Initial values of energy
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

    // Temperature and pressure
    double temperature[length_saved]; // Temperature in K
    double pressure[length_saved];    // Pressure in GPa
    temperature[0] = 0;
    pressure[0] = 0;

    // Equilibration values
    // We found most succes doing the equilibration in
    // 2 steps
    double start_temp = 2;
    double start_pressure = 6;
    double start_temp_2 = 8;
    double start_pressure_2 = 10;
    double timestep_temp = start_temp/dt;
    double timestep_pressure = start_pressure/dt;
    double timestep_temp_2 = start_temp_2/dt;
    double timestep_pressure_2 = start_pressure_2/dt;

    // Values for temp equilibration
    double T_equil = 500.0 + 272.15;           //  K for solid phase
    double tau_t   = 200.0*dt;
    double alpha_t = 1.0;

    // Values for pressure equilibration
    double P_equil      = 1e-4;       // 1e-4 GPa = 1bar 
    double bulk_modulus = 62.0;       // 62 - 102 GPa
    double kappa_p      = 1.0/bulk_modulus;
    double tau_p        = 200.0*dt;
    double alpha_p      = 1.0;

    // Print information before Verlet for convinience
    printf("number of timesteps: %d\n", n_timesteps);
    printf("initial potential energy: %f\n", potential_energy[0]);
    printf("initial kinetic energy: %f\n", kinetic_energy[0]);


    // Time loop
    // Initial calculation of the forces in Verlet algorithm
    get_forces_AL(f, pos, L, n_atoms);

    for(int t = 1; t < n_timesteps + 1; t++){
        
        // Timestepping with Verlet
        verlet_timestep(pos, v, f, n_atoms, dt, m_al, L);

        // Instant values of the system
        // Summing the velocities for kinetic energy
        kinetic_energy[t] = get_kinetic_energy(v, n_atoms, m_al);
        
        // Calculation the potential energy 
        potential_energy[t] = get_energy_AL(pos, L, n_atoms);

        // Calculation of the virial
        virial[t] = get_virial_AL(pos, L, n_atoms);

        // Temperature in K and pressure in GPa
        temperature[t] = (2.0/(3.0*n_atoms))*kinetic_energy[t]/KB; 
        pressure[t] =  (TO_GPA/V[t])*(n_atoms*KB*temperature[t] + virial[t]);
        
        // Rescaling of the velocites for temperature equilibration
        if((t > timestep_temp && t < timestep_pressure) ||
           (t > timestep_temp_2 && t < timestep_pressure_2)){

            alpha_t = 1.0 + (2.0*dt/tau_t)*((T_equil - temperature[t])/temperature[t]);

            for(int i = 0; i < n_atoms; i++){
                for(int j = 0; j < NDIM; j++){
                    v[i][j] = sqrt(alpha_t)*v[i][j];
                }
            }
            // We do not rescale V or L during temp equiilibration
            alpha_p = 1.0;
        }
       
        // Rescaling of the positions for pressure equlibration
        if((t > timestep_pressure && t < timestep_temp_2) ||
           (t > timestep_pressure_2)){

            alpha_p = 1.0 - kappa_p*(dt/tau_p)*(P_equil - pressure[t]);
            
            for(int i = 0; i < n_atoms; i++){
                for(int j = 0; j < NDIM; j++){
                    pos[i][j] = cbrt(alpha_p)*pos[i][j];
                }
            }
            // Rescale size of simulation box
            L = cbrt(alpha_p)*L;
            
            // Update forces the next velocity calculation
            get_forces_AL(f, pos, L, n_atoms);
        }

        // Save the timetrail of the changing volume
        V[t+1] = alpha_p*V[t];
        printf("Saving timestep: %d /%d\n", t, n_timesteps);
    }

    // Printing the needed timetrails
    write_TPV(temperature, pressure, V, time, length_saved, "TPV");
    print_pos(pos, n_atoms, "pos_after_equil");
    print_pos(v, n_atoms, "v_after_equil");

    FILE *ffinal_volume = fopen("final_volume.csv", "w");
    fprintf(ffinal_volume, "%f,", V[length_saved]);
    fclose(ffinal_volume);
}
