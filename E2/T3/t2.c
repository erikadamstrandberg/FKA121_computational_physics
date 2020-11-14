/* 
 * Solveing T2 for E2
 */

// Includes from standard lib
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Our own packages
#include "fft.h"
#include "transform.h"
#include "calc_acc.h"

// Defined variables
#define N_PARTICLES 32
#define PI 3.141592653589

// Prototypes
void write_energy_ord_file(double *U_kin, double *U_pot, 
        double *timesteps, int n_timesteps, int save_every);

void write_energy_norm_file(double E_mode[][N_PARTICLES], int energy_length);

void velocity_verlet(int n_timesteps, int n_particles,
                     double *v, double *q,
                     double *P, double *Q,
                     double *w,
                     double dt,
                     double *m,
		             double kappa, 
                     double alpha,
                     double *U_kin, 
                     double *U_pot,
                     double E_mode[][N_PARTICLES],
                     double trans_matrix[N_PARTICLES][N_PARTICLES],
                     int energy_length,
                     int save_every);


// Main stuff
int main(){
    double Q[N_PARTICLES];
    double P[N_PARTICLES];
    double w[N_PARTICLES];
    double q[N_PARTICLES];
    double v[N_PARTICLES];
    double m[N_PARTICLES];
    double trans_matrix[N_PARTICLES][N_PARTICLES];

    construct_transformation_matrix(trans_matrix, N_PARTICLES);
     
    double kappa = 1.0;
    double alpha = 0.01;
    double E0    = 32.0;

    // Simulation values
    double total_time = 25000;
    double dt         = 0.1;
    
    // How many timesteps to save
    int save_every    = 100;

    int n_timesteps = total_time/dt;
    int energy_length = n_timesteps/save_every;
        
    for (int i = 0; i < N_PARTICLES; i++){
        Q[i] = 0.0;
        P[i] = 0.0;
        m[i] = 1.0;
        w[i] = 2.0*sin((i+1)*PI/(2*(N_PARTICLES+1)));
    }
    P[0] = sqrt(2.0*E0);


    transform_to_normal_modes(trans_matrix, N_PARTICLES, P, v);
    transform_to_normal_modes(trans_matrix, N_PARTICLES, Q, q);

    double timesteps[n_timesteps];
    for (int i = 0; i < n_timesteps; i++) {
        timesteps[i] = i*dt;
    }
 
    double U_kin[energy_length];
    double U_pot[energy_length];
    for (int i = 0; i < energy_length; i++){
        U_kin[i] = 0;
        U_pot[i] = 0;
    }

    double E_mode[energy_length][N_PARTICLES];
    for (int i = 0; i < energy_length; i++) {
        for (int j = 0; j < N_PARTICLES; j++){
            if (i == 0){
                E_mode[0][j] = (1.0/2.0)*(pow(P[j],2)+pow(Q[i]*w[j],2)); 
            } else {
                E_mode[i][j] = 0;
            }
        }
    }

    // Running simulation
    printf("Total number of time steps: %d\n", n_timesteps);
    velocity_verlet(n_timesteps, N_PARTICLES, v, q, P, Q, w, dt, m, kappa, alpha, U_kin, U_pot, E_mode, trans_matrix, energy_length, save_every);


    // After evolving system. 
    write_energy_ord_file(U_kin, U_pot, timesteps, n_timesteps, save_every);

    transform_to_normal_modes(trans_matrix, N_PARTICLES, v, P);
    transform_to_normal_modes(trans_matrix, N_PARTICLES, q, Q);
    

    write_energy_norm_file(E_mode, energy_length);
    return 0;
}

void write_energy_norm_file(double E_mode[][N_PARTICLES], int energy_length){
    FILE *fp = fopen("energy_norm.csv", "w");
    fprintf(fp, "modes\n");

    for(int i = 0; i < energy_length; i++){
        for(int j = 0; j < N_PARTICLES; j++){
            if (j == N_PARTICLES-1){
                fprintf(fp, "%f", E_mode[i][j]);  
            } else {
                fprintf(fp, "%f,", E_mode[i][j]);
            }
        }
        fprintf(fp, "\n");
    }
}


void write_energy_ord_file(double *U_kin, double *U_pot, 
                   double *timesteps, int n_timesteps, int save_every){
    FILE *fp = fopen("energy_ord.csv", "w");
    fprintf(fp, "U_kin,U_pot,time\n");

    int count = 0;
    for(int i = 0; i < n_timesteps; ++i){
        if (i%save_every == 0){
	        fprintf(fp, "%f,%f,%f\n", U_kin[count], U_pot[count], timesteps[i]);
            count += 1;
        }
    }
    fclose(fp);
}

 
void velocity_verlet(int n_timesteps, int n_particles,
                     double *v, double *q,
                     double *P, double *Q,
                     double *w,
                     double dt,
                     double *m,
		             double kappa, 
                     double alpha,
                     double *U_kin, 
                     double *U_pot,
                     double E_mode[][N_PARTICLES],
                     double trans_matrix[N_PARTICLES][N_PARTICLES],
                     int energy_length,
                     int save_every){
    double a[n_particles];
    
    for (int j = 0; j < n_particles; j++) {
            U_kin[0] += m[j]*pow(v[j], 2)/2.0;
    }

    for (int j = 0; j < n_particles+1; j++) {
        if(j == 0) {
            U_pot[0] += pow(q[j], 2)*kappa/2.0;
        } else if(j == n_particles){
            U_pot[0] += pow(q[j-1], 2)*kappa/2.0;
        } else{
            U_pot[0] += pow(q[j]-q[j-1], 2)*kappa/2.0;
        }
    }

    calc_acc(a, q, m, kappa, alpha, n_particles);

    int count = 1;
    for (int i = 1; i < n_timesteps + 1; i++) {

        for (int j = 0; j < n_particles; j++) {
            v[j] += dt * 0.5 * a[j];
        }
        
        /* q(t+dt) */
        for (int j = 0; j < n_particles; j++) {
            q[j] += dt * v[j];
        }
        
        /* a(t+dt) */
        calc_acc(a, q, m, kappa, alpha, n_particles);
        
        /* v(t+dt) */
        for (int j = 0; j < n_particles; j++) {
            v[j] += dt * 0.5 * a[j];
        }
    

        if (i%save_every == 0){
            printf("Saving timestep: %d\n", i);
            /*U_kin(t+dt) */
            for (int j = 0; j < n_particles; j++) {
                U_kin[count] += m[j]*pow(v[j], 2)/2.0;
            }

            /*U_pot(t+dt) */
            for (int j = 0; j < n_particles+1; j++) {
                if(j == 0) {
                    U_pot[count] += pow(q[j], 2)*kappa/2.0 + pow(q[j], 3)*alpha/3.0; 
                } else if(j == n_particles){
                    U_pot[count] += pow(q[j-1], 2)*kappa/2.0 + pow(q[j-1], 3)*alpha/3.0;
                } else{
                    U_pot[count] += pow(q[j]-q[j-1], 2)*kappa/2.0 + pow(q[j]-q[j-1], 3)*alpha/3.0;
                }
            }

            transform_to_normal_modes(trans_matrix, N_PARTICLES, v, P);
            transform_to_normal_modes(trans_matrix, N_PARTICLES, q, Q);
            for (int j = 0; j < N_PARTICLES; j++){ 
                E_mode[count][j] += (1.0/2.0)*(pow(P[j], 2)+pow(Q[j]*w[j], 2)); 
            }
            count += 1;
        }
    }
}

