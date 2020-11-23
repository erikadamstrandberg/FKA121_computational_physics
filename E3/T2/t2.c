// Includes from standard lib
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define PI 3.14159265359

// Main 
int main(){
    int N = 10000;
    double x = 0.0;
    double y = 0.0;
    double gi = 0.0;
    double gi2 = 0.0;
    double I = 0.0;
    double sigma2 = 0.0;

    srand(time(NULL));

    for(int i = 0; i < N; i++){ 
        x = (double) rand()/(double) RAND_MAX;
        y = (1.0/PI)*acos(1.0-2.0*x);
        
        gi += (2.0/PI)*(y*(1-y))/(sin(PI*y));
        gi2 += pow((2.0/PI)*(y*(1-y))/(sin(PI*y)),2);
    }

    I = gi/N;
    sigma2 = gi2/N-pow(gi/N,2);

    printf("%f\n", I);
    printf("%f\n", sigma2);
}
