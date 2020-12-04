// H2b Variational Monto Carlo
// 
// T1 benchmarking

// Standard C libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// GSL for random number generation
#include <gsl/gsl_rng.h>


// Includes from objective files
#include "init_gsl.h"

// Main 
int main(){

    // Initializes GSL random number generation
    gsl_rng* gsl_rand = init_gsl();

    // Pick a random number
    double rand;
    rand = gsl_rng_uniform(gsl_rand);
    printf("%f\n", rand);


    double alpha = 0.1;
    double x1;
    double y1;
    double z1;
    
    double x2;
    double y2;
    double z2;

    double phi_T;
    double energy;

    double r1  = sqrt(pow(x1,2)      + pow(y1,2)      + pow(z1,2));
    double r2  = sqrt(pow(x2,2)      + pow(y2,2)      + pow(z2,2));
    double r12 = sqrt(pow((x2-x1),2) + pow((y2-y1),2) + pow((z2-z1),2));

    // How to calculate equation 6?!
    phi_T = exp(-2*r1)*exp(-2*r2)*(r12/(2*(1+alpha*r12)));
    






    // How to best calculated unit vector
    energy = -4;
    
}


