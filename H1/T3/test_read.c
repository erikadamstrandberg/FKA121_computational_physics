#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define NDIM 3

int main(){

    int nrows = 4;
    double mat[nrows][3];

    FILE *f = fopen("testData.csv", "r");
    int max = 100;
    char buf[max];
    for(int i = 0; i < nrows; i++){
        fgets(buf, max, f);
       // printf("%s", buf);
        //scanf(buf, "%lf,%lf,%lf", &mat[2][0], &mat[2][1], &mat[2][2]);
    }
    fclose(f);

    double x = 0;
    double y = 0;
    double z = 0;

    char *end;
    char str[] = "43.0,23.5,12.0";
    char *temp;
    temp = strtok(str,",");

//    printf("%f", strtod(temp,NULL));

    x = strtod(temp, &end);
    temp = strtok(NULL, ",");
    y = strtod(temp, NULL);
    temp = strtok(NULL, ",");
    z = strtod(temp, NULL);



    //while (temp != NULL){
      //  temp = strtok(NULL, ",");
       // printf("%f", strtod(temp, NULL));
//    }

    // Print mat
    for(int i = 0; i < nrows; i++){
       // printf("%f, %f, %f\n", mat[i][0], mat[i][1], mat[i][2]);
    }

}

void read_initial_position(double *pos, int n_atoms){
     FILE *f = fopen("data/pos_after_equil.csv", "r");
     for(int i = 0; i < n_atoms; i++){
         for(int j = 0; j < NDIM; j++){
            //pos[i][j] = fscanf(f);
            //fprint("atom %d\t %f\n", i, fscanf(f))
         }
     }
}   
