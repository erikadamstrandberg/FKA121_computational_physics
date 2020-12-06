#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


void print_markov(double *E_l, int N,
                  double *x1, double *y1, double *z1, double *x2, double *y2, double *z2,
                  char*filename)
{
    char *filename_csv = malloc((strlen(filename) + 4)*sizeof(char));

    strcpy(filename_csv, filename);
    strcat(filename_csv, ".csv");

    FILE *mc = fopen(filename_csv, "w");
    for(int i = 0; i < N; i++)
    {
       fprintf(mc, "%f,%f,%f,%f,%f,%f,%f\n", x1[i], y1[i], z1[i], x2[i], y2[i], z2[i], E_l[i]);
    }
    fclose(mc);
}
