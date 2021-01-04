#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>

void print_1d_array(double *array, int N, char*filename)
{
    char *filename_csv = malloc((strlen(filename) + 4)*sizeof(char));

    strcpy(filename_csv, filename);
    strcat(filename_csv, ".csv");

    FILE *el = fopen(filename_csv, "w");
    for(int i = 0; i < N; i++)
    {
       fprintf(el, "%f\n", array[i]);
    }
    fclose(el);
}

void print_complex_array(complex *array, int N, char*filename)
{
    char *filename_csv = malloc((strlen(filename) + 4)*sizeof(char));

    strcpy(filename_csv, filename);
    strcat(filename_csv, ".csv");

    FILE *el = fopen(filename_csv, "w");
    for(int i = 0; i < N; i++)
    {
       fprintf(el, "%f, %f\n", creal(array[i]), cimag(array[i]));
    }
    fclose(el);
}

