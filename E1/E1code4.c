/******************************************************************************
 * E1code4
 ******************************************************************************
 * Routine that runs the velocity verlet algorithm
 * Use as template to construct your program!
 */



/*
 * Calculate the acceleration
 * @a - vector that is filled with acceleration
 * @u - vector with the current positions
 * @m - vector with masses
 * @kappa - Spring constant
 * @size_of_u - the size of the position, acceleration and mass array
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_const_mksa.h>
#include <math.h>

#define AMU GSL_CONST_MKSA_UNIFIED_ATOMIC_MASS

void write_qs_file(double *q1, double *q2, double *q3, 
                   double *U_kin, double *U_pot, 
                   double *timesteps, int n_points)
{
    FILE *fp = fopen("timetrail.csv", "w");
    fprintf(fp, "q1,q2,q3,U_kin,U_pot,time\n");
    for(int i = 0; i < n_points; ++i){
	    fprintf(fp, "%f,%f,%f,%f,%f,%f\n", q1[i], q2[i], q3[i], U_kin[i], U_pot[i], timesteps[i]);
    }
    fclose(fp);
}

void calc_acc(double *a, double *u, double *m, double kappa, int size_of_u)
{
    /* Declaration of variables */
    int i;
    
    /* Calculating the acceleration on the boundaries */
    a[0] = kappa*(- u[0] + u[1])/m[0];
    a[size_of_u - 1] = kappa*(u[size_of_u - 2] - u[size_of_u - 1])/m[size_of_u - 1];
    
    /* Calculating the acceleration of the inner points */
    for (i = 1; i < size_of_u - 1; i++){
        a[i] = kappa*(u[i - 1] - 2*u[i] + u[i + 1])/m[i];
    }
}

/*
 * Perform the velocity verlet alogrithm 
 * @n_timesteps - The number of time steps to be performed
 * @n_particles - number of particles in the system
 * @v - array of velocity (Empty allocated array) - sizeof(q_n) = n_timesteps
 * @q_n - position of the n'th atom : sizeof(q_n) = n_timesteps
 * @dt - timestep
 * @m - vector with masses of atoms sizeof(n_particles)
 * @kappa - Spring constant
 */
void velocity_verlet(int n_timesteps, int n_particles, double *v, double *q_1,
		     double *q_2, double *q_3, double dt, double *m,
		     double kappa, 
             double *U_kin, double *U_pot)
{
    double q[n_particles];
    double a[n_particles];
    q[0] = q_1[0];
    q[1] = q_2[0];
    q[2] = q_3[0];
    calc_acc(a, q, m, kappa, n_particles);
    for (int i = 1; i < n_timesteps + 1; i++) {

        if (i%100 == 0) {
            printf("Working on timestep: %d\n", i);
        }

        /* v(t+dt/2) */
        for (int j = 0; j < n_particles; j++) {
            v[j] += dt * 0.5 * a[j];
        }
        
        /* q(t+dt) */
        for (int j = 0; j < n_particles; j++) {
            q[j] += dt * v[j];
        }
        
        /* a(t+dt) */
        calc_acc(a, q, m, kappa, n_particles);
        
        /* v(t+dt) */
        for (int j = 0; j < n_particles; j++) {
            v[j] += dt * 0.5 * a[j];
        }
        
        /*U(t+dt) */
        for (int j = 0; j < n_particles; j++) {
            U_kin[i] += m[j]*v[j]*v[j]/2.0;

            if(j == 0 || j == n_particles-1){
                U_pot[i] += pow(q[j], 2)*kappa/2.0;
            }
            else{
                U_pot[i] += pow(q[j+1]-q[j], 2)*kappa/2.0;
            }
        }

        /* Save the displacement of the three atoms */
        q_1[i] = q[0];
        q_2[i] = q[1];
        q_3[i] = q[2];
    }
}

int main(){

    double total_time = 0.25;
    double dt = 0.0001;
    int n_timesteps = total_time/dt;
    printf("Total number of time steps: %d\n", n_timesteps);
    int n_particles = 3;

    double kappa = 1000e-24/(9649*AMU);

    double v1 = 0.0;
    double v2 = 0.0;
    double v3 = 0.0;
    double v[] = {v1, v2, v3};
    
    double q1[n_timesteps];
    double q2[n_timesteps];
    double q3[n_timesteps];
    q1[0] = 0.01;

    double m1 = 12.0*AMU;
    double m2 = 12.0*AMU;
    double m3 = 12.0*AMU;
    double m[] = {m1, m2, m3};
    
    double U_kin[n_timesteps];
    double U_pot[n_timesteps];
    //U_pot[0] = kappa*pow(q1[0], 2);

    velocity_verlet(n_timesteps, n_particles, v, q1, q2, q3, dt, m, kappa, U_kin, U_pot);

    double timesteps[n_timesteps];
    for (int i = 0; i < n_timesteps; i++) {
        timesteps[i] = i*dt;
    }

    write_qs_file(q1, q2, q3, U_kin, U_pot, timesteps, n_timesteps);
    return 0;
}
