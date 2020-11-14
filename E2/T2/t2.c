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
#define N_PARTICLES 32
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
                   double *timesteps, int n_timesteps, int save_every);


void calc_acc(double *a, double *u, double *m, double kappa, double alpha, int size_of_u);

void velocity_verlet(int n_timesteps, int n_particles,
                     double *v, double *q,
                     double dt,
                     double *m,
		             double kappa, 
                     double alpha,
                     double *U_kin, 
                     double *U_pot,
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
     
    double kappa = 1;
    double alpha = 0;
    double E0 = 32.0;

    // Simulation values
    int save_every = 500;
    double total_time = 25000;
    double dt = 0.1;
    int n_timesteps = total_time/dt;
    
    int energy_length = n_timesteps/save_every;
        
    for (int i = 0; i < N_PARTICLES; i++){
        Q[i] = 0;
        P[i] = 0;
        m[i] = 1;
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

    printf("Total number of time steps: %d\n", n_timesteps);
    velocity_verlet(n_timesteps, N_PARTICLES, v, q, dt, m, kappa, alpha, U_kin, U_pot, save_every);

    write_energy_file(U_kin, U_pot, timesteps, n_timesteps, save_every);

    transform_to_normal_modes(trans_matrix, N_PARTICLES, v, P);
    transform_to_normal_modes(trans_matrix, N_PARTICLES, q, Q);
    
    double E_tot = 0.0;
    double E_mode[N_PARTICLES];
    for (int i = 0; i < N_PARTICLES; i++){
        E_mode[i] = (1.0/2.0)*(pow(P[i],2)+pow(Q[i]*w[i],2));
        E_tot += E_mode[i];
        printf("Mode number %d has energy %.10f\n", i, E_mode[i]);
    }

    printf("Total energy is %.10f\n", E_tot);
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
			       double *from, double *to){
    for(int i = 0; i < n_particles; i++){
	    double sum = 0;
	    for(int j = 0; j < n_particles; j++){
	        sum += from[j] * trans_matrix[i][j];
	    }
	to[i] = sum;
    }
}

void write_energy_file(double *U_kin, double *U_pot, 
                   double *timesteps, int n_timesteps, int save_every){
    FILE *fp = fopen("energy.csv", "w");
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

void calc_acc(double *a, double *u, double *m, double kappa, double alpha, int size_of_u){
    int i;

    /* Calculating the acceleration on the boundaries */
    a[0] = kappa*(-2*u[0] + u[1])/m[0]
          -alpha*(pow(u[0],2)-pow(u[1]-u[0],2))/m[0];
    a[size_of_u - 1] = kappa*(u[size_of_u - 2] - 2*u[size_of_u - 1])/m[size_of_u - 1]
          -alpha*(pow(u[size_of_u-1]-u[size_of_u-2],2)-pow(u[size_of_u-1],2))/m[size_of_u-1];

    /* Calculating the acceleration of the inner points */
    for (i = 1; i < size_of_u - 1; i++){
        a[i] = kappa*(u[i-1]-2*u[i]+u[i+1])/m[i]
              -alpha*(pow(u[i]-u[i-1],2)-pow(u[i+1]-u[i],2))/m[i];
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
            printf("%d", count); 
            /*U_kin(t+dt) */
            for (int j = 0; j < n_particles; j++) {
                U_kin[count] += m[j]*pow(v[j], 2)/2.0;
            }

            /*U_pot(t+dt) */
            for (int j = 0; j < n_particles+1; j++) {
                if(j == 0) {
                    U_pot[count] += pow(q[j], 2)*kappa/2.0;
                } else if(j == n_particles){
                    U_pot[count] += pow(q[j-1], 2)*kappa/2.0;
                } else{
                    U_pot[count] += pow(q[j]-q[j-1], 2)*kappa/2.0;
                }
            }
            count += 1;
        }
    }
}

