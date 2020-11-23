// Includes from standard lib
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void MC(double *I, double *sigma2, int N){
    srand(time(NULL));
    double x;
    double fi = 0.0;
    double fi2 = 0.0;
    
    for(int i = 0; i < N; i++){
        x = (double) rand()/(double) RAND_MAX;
        fi += x*(1.0-x);
        fi2 += pow(x*(1.0-x),2); 
    }
    *I = fi/N;
    *sigma2 = fi2/N-pow(fi/N, 2);
}

void integrate(int N){
    double I = 0.0;
    double sigma2 = 0.0;

    MC(&I, &sigma2, N);
    

    printf("N = %d\n", N);
    printf("I = %f\n", 1.0/6.0);
    printf("I_N = %f\n", I);
    printf("sigma2 = %f\n", sigma2);
    printf("sigma = %f\n", sqrt(sigma2)/sqrt(N));
}

// Main 
int main(){
    int N[] = {10, 100, 1000, 10000};

    for(int i = 0; i < 4; i++){
        integrate(N[i]);
    }
}
