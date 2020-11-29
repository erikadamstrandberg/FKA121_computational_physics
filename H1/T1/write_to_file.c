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

