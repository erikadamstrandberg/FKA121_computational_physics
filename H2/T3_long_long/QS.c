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
    double alpha = 0.1;
    int N_burn   = 10000;
    int N_tot    = 1000000;
    int N_ns     = 10000;
    int N = N_tot - N_ns;
    
    // All values for estimating E_l
    double E_l = 0.0;
    double E_l_mean = 0.0;
    double E_l2_mean = 0.0;
    double E_l_sigma2;
   
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

    // Variables for gradiant decent
    double grad_E;
    double fst_term = 0.0;
    double snd_term = 0.0;
    double grad_ln_phi = 0.0;
    
    // Picking random first step 
    x1 = delta*(gsl_rng_uniform(gsl_rand)-0.5);
    y1 = delta*(gsl_rng_uniform(gsl_rand)-0.5);
    z1 = delta*(gsl_rng_uniform(gsl_rand)-0.5);
    x2 = delta*(gsl_rng_uniform(gsl_rand)-0.5);
    y2 = delta*(gsl_rng_uniform(gsl_rand)-0.5);
    z2 = delta*(gsl_rng_uniform(gsl_rand)-0.5);

    // Initializing the run and finding ns for this alpha
    double ns;
    initialize_mcmc(&ns, &alpha, 
                    &E_l_mean, &E_l2_mean, &fst_term, &snd_term, 
                    &delta, accept, 
                    N_ns, N_burn, 
                    &x1, &y1, &z1, &x2, &y2, &z2, 
                    gsl_rand);
 
    printf("%f\n", ns);
    // Generation of Markov chain with the Metropolis algorithm
    for(int i = 0; i < N; i++)
    {
        metropolis_move(&x1, &y1, &z1, &x2, &y2, &z2, &delta, &alpha, &accept, gsl_rand);
        local_energy(&E_l, &alpha, &x1, &y1, &z1, &x2, &y2, &z2);
        E_l_mean += E_l;
        E_l2_mean += pow(E_l, 2);
        
        grad_alpha_ln_phi(&grad_ln_phi, &alpha, &x1, &y1, &z1, &x2, &y2, &z2); 
        fst_term += E_l*grad_ln_phi;
        snd_term += grad_ln_phi;

    }

    // Calculate estimate of E
    E_l_mean = E_l_mean/N_tot;
    E_l_sigma2 = E_l2_mean/N_tot - pow(E_l_mean, 2);

    // Calculate values for steepest decent
    fst_term = fst_term/N_tot;
    snd_term = E_l_mean*(snd_term/N_tot);
    grad_E = fst_term - snd_term;

    printf("Acceptance ratio: %f\n", (double) accept/((double) N_tot));
    printf("Mean: %f\n", E_l_mean);
    printf("error bar: %f\n", sqrt(ns*E_l_sigma2/N_tot));
}

