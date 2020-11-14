#include <math.h>

#define N_PARTICLES 32
#define PI 3.141592653589

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

