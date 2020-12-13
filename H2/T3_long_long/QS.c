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


void get_ns(double *ns, double *alpha, double *E_l_mean, double *E_l2_mean, double *delta, int accept, int N_ns, int burn_in, 
                 double *x1, double *y1, double *z1, double *x2, double *y2, double *z2, gsl_rng *gsl_rand)
{

    

    double *x1_ns = malloc(N_ns*sizeof(double));
    double *y1_ns = malloc(N_ns*sizeof(double));
    double *z1_ns = malloc(N_ns*sizeof(double));
    double *x2_ns = malloc(N_ns*sizeof(double));
    double *y2_ns = malloc(N_ns*sizeof(double));
    double *z2_ns = malloc(N_ns*sizeof(double));
    double *E_l   = malloc(N_ns*sizeof(double));

    double *mean_E_l  = malloc(sizeof(double));
    double *sigma2_E_l = malloc(sizeof(double));
    
    *mean_E_l = 0.0;
    *sigma2_E_l = 0.0; 

    x1_ns[0] = *x1;
    y1_ns[0] = *y1;
    z1_ns[0] = *z1;
    x2_ns[0] = *x2;
    y2_ns[0] = *y2;
    z2_ns[0] = *z2;

    for(int i = 1; i < burn_in + 1; i++)
    {

        x1_ns[i] = x1_ns[i-1];
        y1_ns[i] = y1_ns[i-1];
        z1_ns[i] = z1_ns[i-1];
        x2_ns[i] = x2_ns[i-1];
        y2_ns[i] = y2_ns[i-1];
        z2_ns[i] = z2_ns[i-1];

        metropolis_move(&x1_ns[i], &y1_ns[i], &z1_ns[i], &x2_ns[i], &y2_ns[i], &z2_ns[i],
                        delta, alpha, &accept, gsl_rand);
    }

    accept = 0;

    x1_ns[0] = x1_ns[burn_in];
    y1_ns[0] = y1_ns[burn_in];
    z1_ns[0] = z1_ns[burn_in];
    x2_ns[0] = x2_ns[burn_in];
    y2_ns[0] = y2_ns[burn_in];
    z2_ns[0] = z2_ns[burn_in];
    
    local_energy(&E_l[0], alpha, &x1_ns[0], &y1_ns[0], &z1_ns[0], &x2_ns[0], &y2_ns[0], &z2_ns[0]);

    for(int i = 1; i < N_ns; i++)
    {

        x1_ns[i] = x1_ns[i-1];
        y1_ns[i] = y1_ns[i-1];
        z1_ns[i] = z1_ns[i-1];
        x2_ns[i] = x2_ns[i-1];
        y2_ns[i] = y2_ns[i-1];
        z2_ns[i] = z2_ns[i-1];

        metropolis_move(&x1_ns[i], &y1_ns[i], &z1_ns[i], &x2_ns[i], &y2_ns[i], &z2_ns[i],
                        delta, alpha, &accept, gsl_rand);
        local_energy(&E_l[i], alpha, &x1_ns[i], &y1_ns[i], &z1_ns[i], &x2_ns[i], &y2_ns[i], &z2_ns[i]);
        *E_l_mean += E_l[i];
        *E_l2_mean += pow(E_l[i], 2);
    }

    mean(mean_E_l, E_l, 0, N_ns);
    sigma2(sigma2_E_l, E_l, 0, N_ns);
    estimate_ns(ns, E_l, *sigma2_E_l, N_ns);
    
    *x1 = x1_ns[N_ns];
    *y1 = y1_ns[N_ns];
    *z1 = z1_ns[N_ns];
    *x2 = x2_ns[N_ns];
    *y2 = y2_ns[N_ns];
    *z2 = z2_ns[N_ns];

    free(x1_ns); x1_ns = NULL;
    free(y1_ns); y1_ns = NULL;
    free(z1_ns); z1_ns = NULL;
    free(x2_ns); x2_ns = NULL;
    free(y2_ns); y2_ns = NULL;
    free(z2_ns); z2_ns = NULL;

}


// Main 
int main()
{
    // Initializes GSL random number generation
    gsl_rng* gsl_rand = init_gsl();

    // Initializing the simulation
    double alpha = 0.25;
    int N_burn = 4000;
    int N_tot = 10000;
    int N_ns = 100000;

    int N = N_tot - N_ns;
    
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

    double ns;
    
    get_ns(&ns, &alpha, &E_l_mean, &E_l2_mean, &delta, accept, N_ns, N_burn, &x1, &y1, &z1, &x2, &y2, &z2, gsl_rand);

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

    E_l_mean = E_l2_mean/N_tot;
    double sigma2 = E_l2_mean/N_tot - pow(E_l_mean, 2);

    fst_term = fst_term/N;
    snd_term = E_l_mean*(snd_term/N);
    grad_E = fst_term - snd_term;

    printf("Acceptance ratio: %f\n", (double) accept/((double) N));
    printf("Mean: %f\n", E_l_mean);
    printf("error bar: %f\n", sqrt(10*sigma2/N));
}
