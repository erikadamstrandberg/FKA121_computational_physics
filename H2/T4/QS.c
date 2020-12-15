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
    int N_burn   = 10000;
    int N_tot    = 60000000;
    int N_ns     = 10000000;
    int N = N_tot - N_ns;
    int N_opt = 10;

    // Setting MC step length  
    double delta = 1.0;
    
    // All values for estimating E_l
    double E_l = 0.0;
    double E_l_mean = 0.0;
    double E_l2_mean = 0.0;
    double E_l_sigma2 = 0.0;
   
    double ns = 0.0;

    // Variables for generating the Markovs chain
    double x1;
    double y1;
    double z1;    
    double x2;
    double y2;
    double z2;

    // Recordning the acceptance ratio
    int accept = 0;

    // Variables for gradiant decent
    double grad_E = 0.0;
    double fst_term = 0.0;
    double snd_term = 0.0;
    double grad_ln_phi = 0.0;

    // Array for storing all the optimized alphas
    double alpha[N_opt];

    // Setting a starting point for optimization
    // and beta for damped steepest decent
    alpha[0] = 0.14;
    double beta = 1.0;

    // Looping over how many steps that we want to optimize
    for(int p = 1; p < N_opt + 1; p++)
    { 

        // Resetting all variables from last optimization
        E_l = 0.0;
        E_l_mean = 0.0;
        E_l2_mean = 0.0;
        E_l_sigma2 = 0.0;

        ns = 0.0;
        accept = 0;
        
        grad_E = 0.0;
        fst_term = 0.0;
        snd_term = 0.0;
        grad_ln_phi = 0.0;
    
        // Picking random first step 
        x1 = delta*(gsl_rng_uniform(gsl_rand)-0.5);
        y1 = delta*(gsl_rng_uniform(gsl_rand)-0.5);
        z1 = delta*(gsl_rng_uniform(gsl_rand)-0.5);
        x2 = delta*(gsl_rng_uniform(gsl_rand)-0.5);
        y2 = delta*(gsl_rng_uniform(gsl_rand)-0.5);
        z2 = delta*(gsl_rng_uniform(gsl_rand)-0.5);

        // Initializing the run and finding ns for this alpha
        initialize_mcmc(&ns, &alpha[p-1], 
                        &E_l_mean, &E_l2_mean, &fst_term, &snd_term, 
                        &delta, &accept, 
                        N_ns, N_burn, 
                        &x1, &y1, &z1, &x2, &y2, &z2, 
                        gsl_rand);

        // Generation of Markov chain with the Metropolis algorithm
        for(int i = 0; i < N; i++)
        {
            metropolis_move(&x1, &y1, &z1, &x2, &y2, &z2, &delta, &alpha[p-1], &accept, gsl_rand);
            local_energy(&E_l, &alpha[p-1], &x1, &y1, &z1, &x2, &y2, &z2);
            E_l_mean += E_l;
            E_l2_mean += pow(E_l, 2);
            
            grad_alpha_ln_phi(&grad_ln_phi, &alpha[i-1], &x1, &y1, &z1, &x2, &y2, &z2); 
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
        
        alpha[p] = alpha[p-1] - grad_E/pow(p, beta);

        printf("Acceptance ratio: %f\n", (double) accept/((double) N_tot));
        printf("Mean: %f\n", E_l_mean);
        printf("error bar: %f\n", sqrt(ns*E_l_sigma2/N_tot));
        printf("ns: %f\n", ns);
        printf("\nalpha: %f\n\n", alpha[p]);
    }
}

