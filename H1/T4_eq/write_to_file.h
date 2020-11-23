#define NDIM 3

void print_pos(double pos[][NDIM], int n_atoms, char *filename);
void write_TPV(double *temperature, double *pressure, double *V, double *time, int length_saved, char *filename);
