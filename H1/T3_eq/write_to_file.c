#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define NDIM 3

void print_pos(double pos[][NDIM], int n_atoms, char *filename){
    char *filename_csv = malloc((strlen(filename) + 4)*sizeof(char));

    strcpy(filename_csv, filename);
    strcat(filename_csv, ".csv");

    FILE *fpos = fopen(filename_csv, "w");
    for(int i = 0; i < n_atoms; i++){
        for(int j = 0; j < NDIM; j++){
            if(j==2){
                fprintf(fpos, "%f\n", pos[i][j]);
            } else {
                fprintf(fpos, "%f,", pos[i][j]);
            }
        }
    }
    fclose(fpos);
    free(filename_csv);
}

void write_TPV(double *temperature, double *pressure, double *V, double *time, int length_saved, char *filename){
    int length_filename = strlen(filename);
    char *filename_csv = malloc((length_filename + 4)*sizeof(char));
    strcpy(filename_csv, filename);
    strcat(filename_csv, ".csv");

    FILE *ftemp_press = fopen(filename_csv, "w");
    fprintf(ftemp_press, "temp[K],pressure[Pa]\n");
    for (int i = 0; i < length_saved; i++){
        fprintf(ftemp_press, "%f,%f,%f,%f\n", temperature[i], pressure[i], V[i], time[i]);
    }
    fclose(ftemp_press);
    free(filename_csv);
}

void write_energy(double *kinetic_energy, double *potential_energy, double *time, int length_saved){
    FILE *fp = fopen("energy.csv", "w");
    fprintf(fp, "kinetic_energy,potential_energy,time\n");
    for(int t = 0; t < length_saved; t++){
        fprintf(fp, "%f,%f,%f\n", kinetic_energy[t], potential_energy[t], time[t]);
    }
}
