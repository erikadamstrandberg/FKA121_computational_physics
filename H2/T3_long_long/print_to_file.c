#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


void print_markov(int N,
                  double *x1, double *y1, double *z1, double *x2, double *y2, double *z2,
                  char*filename)
{
    char *filename_csv = malloc((strlen(filename) + 4)*sizeof(char));

    strcpy(filename_csv, filename);
    strcat(filename_csv, ".csv");

    FILE *mc = fopen(filename_csv, "w");
    for(int i = 0; i < N; i++)
    {
       fprintf(mc, "%f,%f,%f,%f,%f,%f\n", x1[i], y1[i], z1[i], x2[i], y2[i], z2[i]);
    }
    fclose(mc);
}

void print_1d_array(double *E_l, int N,
                  char*filename)
{
    char *filename_csv = malloc((strlen(filename) + 4)*sizeof(char));

    strcpy(filename_csv, filename);
    strcat(filename_csv, ".csv");

    FILE *el = fopen(filename_csv, "w");
    for(int i = 0; i < N; i++)
    {
       fprintf(el, "%f\n", E_l[i]);
    }
    fclose(el);
}

void print_current_state(double *x1, double *y1, double *z1, 
                         double *x2, double *y2, double *z2,
                         char*filename)
{
    char *filename_csv = malloc((strlen(filename) + 4)*sizeof(char));

    strcpy(filename_csv, filename);
    strcat(filename_csv, ".csv");

    FILE *cs = fopen(filename_csv, "w");
    fprintf(cs, "%f,%f,%f,%f,%f,%f\n", *x1, *y1, *z1, *x2, *y2, *z2);
    fclose(cs);
}
