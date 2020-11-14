#define N_PARTICLES 32

void construct_transformation_matrix(
                            double trans_matrix[N_PARTICLES][N_PARTICLES], 
                            int n_particles);
void transform_to_normal_modes(
                            double trans_matrix[N_PARTICLES][N_PARTICLES],
			                int n_particles,
			                double *q, double *Q);
