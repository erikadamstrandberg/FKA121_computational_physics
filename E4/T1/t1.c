// Main function to solve E4
//
// Task 1. Implementing the BD3 algorithm
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_randist.h>


// Main 
int main(){
    
    const gsl_rng_type *T;
    gsl_rng *Q;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    Q = gsl_rng_alloc(T);
    gsl_rng_set(Q, time(NULL));
    
    int N = 10000;
    double x[N];
    for (int i = 0; i < N; i++){
        x[i] = gsl_ran_gaussian(Q, 1);
    }
 
    FILE *fp = fopen("gaussian.csv", "w");
    for(int i = 0; i < N; i++){
        fprintf(fp, "%f\n", x[i]);
    }
    fclose(fp);

}
