// Includes from standard lib
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define PI 3.14159265359

void MC(int N){

    double x[N];
    double y[N];
    double gi = 0.0;
    double gi2 = 0.0;
    double I = 0.0;
    double sigma2 = 0.0;

    srand(time(NULL));

    for(int i = 0; i < N; i++){ 
        x[i] = (double) rand()/(double) RAND_MAX;
        y[i] = (1.0/PI)*acos(1.0-2.0*x[i]);
        
        gi += (2.0/PI)*(y[i]*(1-y[i]))/(sin(PI*y[i]));
        gi2 += pow((2.0/PI)*(y[i]*(1-y[i]))/(sin(PI*y[i])),2);
    }

    I = gi/N;
    sigma2 = gi2/N-pow(gi/N,2);

    printf("N = %d\n", N);
    printf("I = %f\n", 1.0/6.0);
    printf("I_n = %f\n", I);
    printf("sigma2 = %f\n", sigma2);
    printf("sigma = %f\n\n", sqrt(sigma2)/sqrt(N));
    printf("I = %f +- %f\n\n", I, sqrt(sigma2)/sqrt(N));

    char filename[100];
    char name[] = "dist_";
    char number[100];
    char csv[] = ".csv";
    sprintf(number, "%d", N);

    strcpy(filename, name);
    strcat(filename, number);
    strcat(filename, csv);

    FILE *fp = fopen(filename, "w");
    fprintf(fp, "x,y\n");
    for(int i = 0; i < N; i++){
        fprintf(fp, "%f,%f\n", x[i], y[i]);
    }
    fclose(fp);

}


// Main 
int main(){
    int N[] = {10,100,1000,10000,100000};
    int number_of_N = 5;
    for(int i = 0; i < number_of_N; i++){
        MC(N[i]);
    }
}
