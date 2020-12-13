// H2b Variational Monto Carlo
// 
// T1 benchmarking

// Standard C libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// GSL for random number generation
#include <gsl/gsl_rng.h>

// Includes from objective files
#include "init_gsl.h"
#include "mcmc_sampling.h"
#include "print_to_file.h"
#include "error_estimate.h"


// Main 
int main()
{
    // Initializes GSL random number generation
    gsl_rng* gsl_rand = init_gsl();

    // Initializing the simulation
    double alpha = 0.05;
    int N_burn = 4000;
    int N = 120000;
    double E_l[N];

    // Variables for generating the Markovs chain
    double x1;
    double y1;
    double z1;
    
    double x2;
    double y2;
    double z2;

    // Setting MC step length  
    double delta = 1.0;

    // Recordning the acceptance ratio
    int accept = 0;

    // Picking random first step 
    x1 = delta*(gsl_rng_uniform(gsl_rand)-0.5);
    y1 = delta*(gsl_rng_uniform(gsl_rand)-0.5);
    z1 = delta*(gsl_rng_uniform(gsl_rand)-0.5);
    x2 = delta*(gsl_rng_uniform(gsl_rand)-0.5);
    y2 = delta*(gsl_rng_uniform(gsl_rand)-0.5);
    z2 = delta*(gsl_rng_uniform(gsl_rand)-0.5);

    // run burn in
    for(int i = 0; i < N_burn; i++)
    {
        metropolis_move(&x1, &y1, &z1, &x2, &y2, &z2, &delta, &alpha,  &accept, gsl_rand);
    }
    // Generation of Markov chain with the Metropolis algorithm
    for(int i = 0; i < N; i++)
    {
        metropolis_move(&x1, &y1, &z1, &x2, &y2, &z2, &delta, &alpha,  &accept, gsl_rand);
        local_energy(&E_l[i], &alpha, &x1, &y1, &z1, &x2, &y2, &z2);
    }

    printf("Acceptance ratio: %f\n", (double) accept/((double) N));
    

    // writes current state and the local energy to file
    print_current_state(&x1, &y1, &z1, &x2, &y2, &z2, "current_state");
    print_1d_array(E_l, N, "local_energy");
    
    // Calculation of mean and variance of E_l
    double mean_E_l = 0.0;
    double sigma2_E_l = 0.0;
    mean(&mean_E_l, E_l, 0, N);
    sigma2(&sigma2_E_l, E_l, 0, N);


    double ns;
    estimate_ns(&ns, E_l, sigma2_E_l, N);
    

    // Prints the error bars for the calculations
    FILE *estimates = fopen("finalvalue.csv", "w"); 
    fprintf(estimates, "%f,%f,%f,%d,%f\n", mean_E_l, sigma2_E_l, alpha, N, ns);
    fclose(estimates);

    printf("Mean of E_l: %f\n", mean_E_l);
    printf("Variance of E_l: %f\n", sigma2_E_l);
    printf("Statistical inefficiency: %f\n\n", ns);
    printf("Estimate: %f +- %f\n",  mean_E_l, sqrt(sigma2_E_l)/(sqrt(N/ns)));
}

