// H2b Variational Monto Carlo
// 
// The final source code for H2b

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
    int N_tot    = 50000000;
    int N_ns     = 100000;
    int N = N_tot - N_ns;
    int N_opt = 4;
    int N_run_per_alpha = 200;

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

    // Estimates for each alpha run
    double E_est = 0.0;
    double E_est2 = 0.0;
    double fst_term_est = 0.0;
    double snd_term_est = 0.0;

    // Array for storing all the optimized alphas
    double alpha[N_opt];

    // Setting a starting point for optimization
    // and beta for damped steepest decent
    alpha[0] = 0.143;
    double beta = 1.0;

    // Opening a file to save values during the run
    char filename_csv[80];
    char alpha_buffer[10];
    sprintf(alpha_buffer, "%.3f", alpha[0]);
    strcpy(filename_csv, "data/steepest_decent_fine_");
    strcat(filename_csv, alpha_buffer);
    strcat(filename_csv, ".csv");
    FILE *fp = fopen(filename_csv, "w");

    // Looping over how many steps that we want to optimize
    for(int p = 1; p < N_opt + 1; p++)
    { 
        // Resetting values for one value of alpha
        E_est = 0.0;
        E_est2 = 0.0;
        fst_term_est = 0.0;
        snd_term_est = 0.0;

        // Running the independet runs for one value of alpha
        for(int a = 0; a < N_run_per_alpha; a++)
        {
            // Resetting all variables from last run
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
                
                // Storing value for mean and var of the energy 
                E_l_mean  += E_l;
                E_l2_mean += pow(E_l, 2);
                
                // Storing values for estimate of the gradient
                grad_alpha_ln_phi(&grad_ln_phi, &alpha[p-1], &x1, &y1, &z1, &x2, &y2, &z2); 
                fst_term += E_l*grad_ln_phi;
                snd_term += grad_ln_phi;

            }
            // Estimates for the current alpha[p-1]
            // Calculate estimate of E
            E_l_mean = E_l_mean/N_tot;
            E_l_sigma2 = E_l2_mean/N_tot - pow(E_l_mean, 2);

            // Calculate estimates for the gradient
            fst_term = fst_term/N_tot;
            snd_term = snd_term/N_tot;

            // Estimate for entire run
            E_est  += E_l_mean;
            E_est2 += pow(E_l_mean, 2);
            fst_term_est += fst_term;
            snd_term_est += snd_term;
            
            // Printing to know that the run i alive
            printf("Acceptance ratio: %f\n", (double) accept/((double) N_tot));
            printf("Mean: %f\n", E_l_mean);
            printf("error bar: %f\n", sqrt(ns*E_l_sigma2/N_tot));
            printf("ns: %f\n", ns);
            printf("\nalpha: %f\n\n", alpha[p-1]);
        }
        
        // Calculating all the estimates all the independent
        // run for this alpha
        E_est = E_est/N_run_per_alpha;
        E_est2 = E_est2/N_run_per_alpha;
        fst_term_est = fst_term_est/N_run_per_alpha;
        snd_term_est = snd_term_est/N_run_per_alpha;
       
        // The gradient to iterate to the next alpha
        grad_E = 2.0*(fst_term_est - E_est*snd_term_est);

        // Next alpha to calculate
        alpha[p] = alpha[p-1] - grad_E/pow(p, beta);

        // Printing all values
        fprintf(fp, "%.14f,%.14f,%d,%f,%f\n", E_est, E_est2, N_run_per_alpha, grad_E, alpha[p-1]);
    }
    fclose(fp);
}

