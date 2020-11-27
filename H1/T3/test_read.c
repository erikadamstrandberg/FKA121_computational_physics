#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define NDIM 3
#define N_ATOMS 256

int main(){

    char *temp;
    double pos[N_ATOMS][NDIM];


    FILE *fp = fopen("../T3_eq/data/pos_500K_1fs_1.csv", "r");
    char buf_value[10000];
    char *end;
    int skip_header = 0;

    int atom_read = 0;
    while(fgets(buf_value, N_ATOMS + 1, fp)){
        if (skip_header != 0){
            temp = strtok(buf_value,",");
            pos[atom_read][0] = strtod(temp, &end);
            temp = strtok(NULL, ",");
            pos[atom_read][1] = strtod(temp, NULL);
            temp = strtok(NULL, ",");
            pos[atom_read][2] = strtod(temp, NULL);
            atom_read += 1;
        }
        skip_header += 1;
    }
    fclose(fp);

    for (int i = 0; i < N_ATOMS; i++){
        printf("%f,%f,%f\n", pos[i][0], pos[i][1], pos[i][2]);
    }
}
    



