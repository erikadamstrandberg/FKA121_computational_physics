/* 
 * Solveing T2 for E2
 */

// Includes from standard lib
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_const_mksa.h>


// Our own packages
#include "fft.h"


// Defined variables
#define N_PARTICLES 3
#define PI 3.141592653589
#define AMU GSL_CONST_MKSA_UNIFIED_ATOMIC_MASS


// Prototypes
void construct_transformation_matrix(
                            double trans_matrix[N_PARTICLES][N_PARTICLES], 
                            int n_particles);

void transform_to_normal_modes(
                            double trans_matrix[N_PARTICLES][N_PARTICLES],
			                int n_particles,
			                double *q, double *Q);

void write_energy_file(double *U_kin, double *U_pot, 
                   double *timesteps, int n_timesteps){
    FILE *fp = fopen("energy.csv", "w");
    fprintf(fp, "U_kin,U_pot,time\n");
    for(int i = 0; i < n_timesteps; ++i){
	    fprintf(fp, "%f,%f,%f\n", U_kin[i], U_pot[i], timesteps[i]);
    }
    fclose(fp);
}

void calc_acc(double *a, double *u, double *m, double kappa, double alpha, int size_of_u){
    int i;
    
    /* Calculating the acceleration on the boundaries */
    a[0] = kappa*(-2*u[0] + u[1])/m[0];
    a[size_of_u - 1] = kappa*(u[size_of_u - 2] - 2*u[size_of_u - 1])/m[size_of_u - 1];

    /* Calculating the acceleration of the inner points */
    for (i = 1; i < size_of_u - 1; i++){
        a[i] = kappa*(u[i - 1] - 2*u[i] + u[i + 1])/m[i];
    }
}

void velocity_verlet(int n_timesteps, int n_particles,
                     double *v, double *q,
                     double dt,
                     double *m,
		             double kappa, 
                     double alpha,
                     double *U_kin, 
                     double *U_pot,
                     int save_every){

    double a[n_particles];
    
    for (int j = 0; j < n_particles; j++) {
            U_kin[0] += m[j]*pow(v[j], 2)/2.0;
    }

    /*U_pot(t+dt) */
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
                U_kin[i] += m[j]*pow(v[j], 2)/2.0;
            }

            /*U_pot(t+dt) */
            for (int j = 0; j < n_particles+1; j++) {
                if(j == 0) {
                    U_pot[i] += pow(q[j], 2)*kappa/2.0;
                } else if(j == n_particles){
                    U_pot[i] += pow(q[j-1], 2)*kappa/2.0;
                } else{
                    U_pot[i] += pow(q[j]-q[j-1], 2)*kappa/2.0;
                }
            }
        }
    }
}
// Main stuff
int main(){
    double trans_matrix[N_PARTICLES][N_PARTICLES];
    double q[N_PARTICLES];
    double Q[N_PARTICLES];

    double v[N_PARTICLES];
    double m[N_PARTICLES];

    for (int i = 0; i < N_PARTICLES; i++){
        q[i] = 0;
        m[i] = 12.0/9649.0;
        v[i] = 0;
    }

    q[0] = 0.01;
    q[1] = 0.005;
    q[2] = -0.005;

    int save_every = 1;

    construct_transformation_matrix(trans_matrix, N_PARTICLES);

    double total_time = 1;
    double dt = 0.0001;
    int n_timesteps = total_time/dt;
    printf("Total number of time steps: %d\n", n_timesteps);

    double kappa = 1600.0*1e-24/(9649.0*AMU);
    double alpha = 0;
    printf("kappa = %f\n", kappa);

    // Time scaling comes from m and kappa 
    // t = sqrt(mk)
    //
    // We set the length scale to ångstrom
    // This gives the scaling factor for energy as
    
    double U_kin[n_timesteps/save_every];
    double U_pot[n_timesteps/save_every];

    for (int i = 0; i < n_timesteps; i++){
        U_kin[i] = 0;
        U_pot[i] = 0;
    }

    velocity_verlet(n_timesteps, N_PARTICLES, v, q, dt, m, kappa, alpha, U_kin, U_pot, save_every);

    double timesteps[n_timesteps];
    for (int i = 0; i < n_timesteps; i++) {
        timesteps[i] = i*dt;
    }
 
    write_energy_file(U_kin, U_pot, timesteps, n_timesteps);
    transform_to_normal_modes(trans_matrix, N_PARTICLES, q, Q);
    return 0;
}


// Functions
/*
 * trans_matrix[N_PARTICLES][N_PARTICLES]: empty allocated array which
 * will be filled with sine transformation matrix
 * N_PARTICLES: number of particles in system
 */
void construct_transformation_matrix(
    double trans_matrix[N_PARTICLES][N_PARTICLES], int n_particles){
    double factor = 1 / ((double)n_particles + 1);
    for(int i = 0; i < n_particles; i++){
	    for(int j = 0; j < n_particles; j++){
	        trans_matrix[i][j] = sqrt(2 * factor)
			    	 * sin((j + 1) * (i + 1) * PI * factor);
	    }
    }
}

/*
 * Transformation matrix constucted in above function
 * q cartesian coordinate of paricles
 * Q output normal modes coordinate
 * N_PARTICLES is number of particles in system
 */
void transform_to_normal_modes(double trans_matrix[N_PARTICLES][N_PARTICLES],
			       int n_particles,
			       double *q, double *Q){
    for(int i = 0; i < n_particles; i++){
	    double sum = 0;
	    for(int j = 0; j < n_particles; j++){
	        sum += q[j] * trans_matrix[i][j];
	    }
	Q[i] = sum;
    }
}

