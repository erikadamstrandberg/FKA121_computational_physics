#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define NDIM 3

void print_pos(double pos[][NDIM], int n_atoms, char *filename){
    int length_filename = strlen(filename);
    char *filename_csv = malloc((length_filename + 4)*sizeof(char));

    strcpy(filename_csv, filename);
    strcat(filename_csv, ".csv");

    FILE *fpos = fopen(filename_csv, "w");
    fprintf(fpos, "x,y,z\n");
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

void write_temp_pressure(double temperature, double pressure, char *filename){
    int length_filename = strlen(filename);
    char *filename_csv = malloc((length_filename + 4)*sizeof(char));
    strcpy(filename_csv, filename);
    strcat(filename_csv, ".csv");

    FILE *ftemp_press = fopen(filename_csv, "w");
    fprintf(ftemp_press, "temp[K],pressure[Pa]\n");
    fprintf(ftemp_press, "%f,%f", temperature, pressure);
    fclose(ftemp_press);
    free(filename_csv);
}

