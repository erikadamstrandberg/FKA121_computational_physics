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

void grad_alpha_ln_phi(double *grad_ln_phi, double *alpha, 
                       double *x1, double *y1, double *z1, double *x2, double *y2, double *z2)
{
    // Allocating variables
    double *r12 = malloc(sizeof(double));

    *r12 = sqrt(pow((*x2-*x1),2) + pow((*y2-*y1),2) + pow((*z2-*z1),2));
    *grad_ln_phi = -pow(*r12,2)/(2.0*pow(1.0+(*alpha)*(*r12),2));
    
    // Free all variables
    free(r12);
    r12 = NULL;
}

// Main 
int main()
{
    // Initializes GSL random number generation
    gsl_rng* gsl_rand = init_gsl();

    // Initializing the simulation
    double alpha = 0.25;
    int N_burn = 1000;
    int N = 10000;
    double E_l_mean = 0.0;
    double E_l2_mean = 0.0;

    double E_l = 0.0;
    

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

    // Variables for gradiant decent
    double grad_E;
    double fst_term = 0.0;
    double snd_term = 0.0;
    double grad_ln_phi = 0.0;
    
    // run burn in
    for(int i = 0; i < N_burn; i++)
    {
        metropolis_move(&x1, &y1, &z1, &x2, &y2, &z2, &delta, &alpha,  &accept, gsl_rand);
    }

    accept = 0.0;

    // Generation of Markov chain with the Metropolis algorithm
    for(int i = 0; i < N; i++)
    {
        metropolis_move(&x1, &y1, &z1, &x2, &y2, &z2, &delta, &alpha,  &accept, gsl_rand);
        local_energy(&E_l, &alpha, &x1, &y1, &z1, &x2, &y2, &z2);
        E_l_mean += E_l;
        E_l2_mean += pow(E_l, 2);
        
        grad_alpha_ln_phi(&grad_ln_phi, &alpha, &x1, &y1, &z1, &x2, &y2, &z2); 
        fst_term += E_l*grad_ln_phi;
        snd_term += grad_ln_phi;

    }

    E_l_mean = E_l_mean/N;
    double sigma2 = E_l2_mean/N - pow(E_l_mean, 2);

    fst_term = fst_term/N;
    snd_term = E_l_mean*(snd_term/N);
    grad_E = fst_term - snd_term;

    printf("Acceptance ratio: %f\n", (double) accept/((double) N));
    printf("Mean: %f\n", E_l_mean);
    printf("error bar: %f\n", sqrt(10*sigma2/N));
}
