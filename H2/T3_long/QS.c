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
#include "local_energy.h"
#include "print_to_file.h"
#include "error_estimate.h"

// Main 
int main()
{
    // Initializes GSL random number generation
    gsl_rng* gsl_rand = init_gsl();

    // Initializing the simulation
    double alpha = 0.05;
    int N_tot = 120000;
    int burn_in = 4000;
    int N = N_tot - burn_in;
    double E_l[N];

    // Variables for generating the Markovs chain
    double x1;
    double y1;
    double z1;
    double x1_t;
    double y1_t;
    double z1_t;
    
    double x2;
    double y2;
    double z2;
    double x2_t;
    double y2_t;
    double z2_t;

    // Variables for storing the weight function
    double w;
    double w_t;

    // Setting the step length  
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

    // Generation of Markov chain with the Metropolis algorithm
    for(int i = 0; i < N; i++)
    {
        x1_t = x1 + delta*(gsl_rng_uniform(gsl_rand)-0.5);
        y1_t = y1 + delta*(gsl_rng_uniform(gsl_rand)-0.5);
        z1_t = z1 + delta*(gsl_rng_uniform(gsl_rand)-0.5);
        x2_t = x2 + delta*(gsl_rng_uniform(gsl_rand)-0.5);
        y2_t = y2 + delta*(gsl_rng_uniform(gsl_rand)-0.5);
        z2_t = z2 + delta*(gsl_rng_uniform(gsl_rand)-0.5);

        weight(&w, &alpha, &x1, &y1, &z1, &x2, &y2, &z2);  
        weight(&w_t, &alpha, &x1_t, &y1_t, &z1_t, &x2_t, &y2_t, &z2_t);

        if(w_t > w || w_t/w > gsl_rng_uniform(gsl_rand))
        {
           x1 = x1_t; 
           y1 = y1_t; 
           z1 = z1_t; 
           x2 = x2_t; 
           y2 = y2_t; 
           z2 = z2_t;
           accept += 1;
        } 
        
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

