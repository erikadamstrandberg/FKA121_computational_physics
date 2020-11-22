#define NDIM 3

void print_pos(double pos[][NDIM], int n_atoms, char *filename);
void write_temp_pressure(double *temperature, double *pressure, double *time, int length_saved, char *filename);
