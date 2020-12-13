#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include "mcmc_sampling.h"

void mean(double *mean, double *array, int start, int stop)
{
    *mean = 0.0;
    for(int i = start; i < stop; i++)
    {
        *mean += array[i];
    }
    *mean = *mean/(stop-start);
}

void sigma2(double *sigma2, double *array, int start, int stop)
{
    double *mean_array = malloc(sizeof(double));
    mean(mean_array, array, start, stop);

    *sigma2 = 0.0;
    for(int i = start; i < stop; i++)
    {
        *sigma2 += pow(array[i], 2);
    }
    *sigma2 = *sigma2/(stop-start) - pow(*mean_array,2);
    free(mean_array); mean_array = NULL;
}

void correlation_function(double *phi, double *E_l, double *sigma2_E_l, int N, int k)
{
    double *mean_forward = malloc(sizeof(double));
    double *mean_lag     = malloc(sizeof(double));

    for(int i = 0; i < k; i++)
    {
        mean(mean_forward, E_l, i, N);
        mean(mean_lag, E_l, 0, N-i);
        for(int j = 0; j < N-i; j++)
        {
            phi[i] += (E_l[i+j] - (*mean_forward))*(E_l[j] - (*mean_lag));
        }
        phi[i] = phi[i]/(*sigma2_E_l*(N-i));
    }

    free(mean_forward); mean_forward = NULL;
    free(mean_lag); mean_lag = NULL;
}

void block_averaging(double *ns, double *E_l, double *sigma2_E_l, int N, int B)
{ 
    double *Mb = malloc(sizeof(double));
    double *sigma2_each_block = malloc(sizeof(double));
    for(int i = 1; i <= B; i++)
    {
        *Mb = N/i;
        double *F = malloc((*Mb)*sizeof(double));
        for(int j = 0; j < *Mb; j++) 
        {
            F[j] = 0.0;
            for(int k = 0; k < i; k++)
            {
                F[j] += E_l[k+j*i]; 
            }
            F[j] = F[j]/i;
        }
        sigma2(sigma2_each_block, F, 0, *Mb);

        ns[i-1] = i*(*sigma2_each_block)/(*sigma2_E_l);
        free(F);
        F = NULL;
    }
    free(Mb);
    free(sigma2_each_block);
    Mb = NULL;
    sigma2_each_block = NULL;
}

void estimate_ns(double *ns, double *E_l, double sigma2_E_l, int N)
{

    // Caculate and writes the correlation function
    int *k = malloc(sizeof(int));
    *k = 30;
    double *phi = malloc((*k)*sizeof(double));    
    correlation_function(phi, E_l, &sigma2_E_l, N, *k);
    
    // Calculates and writes the block averaging
    int *B = malloc(sizeof(int));
    *B = 800;
    double *ns_b = malloc((*B)*sizeof(double));
    block_averaging(ns_b, E_l, &sigma2_E_l, N, *B);
    
    // Estimate ns from correlation function
    double *limit = malloc(sizeof(double));
    *limit = 0.135;
    double *ns_from_corr = malloc(sizeof(double));
    
    int *index_after = malloc(sizeof(int));
    *index_after = 0;
    while (phi[(*index_after)] > *limit)
    {
        *index_after += 1;
    }
    
    *ns_from_corr = (*index_after)*(phi[(*index_after) - 1]- (*limit))/(phi[(*index_after) - 1] - phi[(*index_after)]) + 
                  ((*index_after)-1)*((*limit)-phi[(*index_after)])/(phi[(*index_after) - 1]-phi[(*index_after)]); 

    // Estimate from block average
    double *ns_from_block = malloc(sizeof(double));
    *ns_from_block = 0.0;
    int *blocks_from_end = malloc(sizeof(int));
    *blocks_from_end = 200;
    for(int i = (*B) - (*blocks_from_end); i < (*B); i++)
    {
        *ns_from_block += ns_b[i];
    }
    *ns_from_block = (*ns_from_block)/(*blocks_from_end);

    if (*ns_from_block > *ns_from_corr)
    {
        *ns = *ns_from_block;
    }
    else 
    {
        *ns = *ns_from_corr;
    }

    free(k); k = NULL;
    free(phi); phi = NULL;
    free(B); B = NULL;
    free(ns_b); ns_b = NULL;
    free(limit); limit = NULL;
    free(ns_from_corr); ns_from_corr = NULL;
    free(index_after);  index_after = NULL;
    free(ns_from_block); ns_from_block = NULL;
    free(blocks_from_end); blocks_from_end = NULL;
}

void initialize_mcmc(double *ns, double *alpha,
            double *E_l_mean, double *E_l2_mean,
            double *fst_term, double *snd_term,
            double *delta, int accept, int N_ns, int burn_in, 
            double *x1, double *y1, double *z1, double *x2, double *y2, double *z2,
            gsl_rng *gsl_rand)
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
    double *grad_ln_phi = malloc(sizeof(double));
    
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
        
        grad_alpha_ln_phi(grad_ln_phi, alpha, &x1_ns[i], &y1_ns[i], &z1_ns[i], &x2_ns[i], &y2_ns[i], &z2_ns[i]); 
        *fst_term += E_l[i]*(*grad_ln_phi);
        *snd_term += *grad_ln_phi;

    }

    mean(mean_E_l, E_l, 0, N_ns);
    sigma2(sigma2_E_l, E_l, 0, N_ns);
    estimate_ns(ns, E_l, *sigma2_E_l, N_ns);
    
    *x1 = x1_ns[N_ns - 1];
    *y1 = y1_ns[N_ns - 1];
    *z1 = z1_ns[N_ns - 1];
    *x2 = x2_ns[N_ns - 1];
    *y2 = y2_ns[N_ns - 1];
    *z2 = z2_ns[N_ns - 1];

    free(x1_ns); x1_ns = NULL;
    free(y1_ns); y1_ns = NULL;
    free(z1_ns); z1_ns = NULL;
    free(x2_ns); x2_ns = NULL;
    free(y2_ns); y2_ns = NULL;
    free(z2_ns); z2_ns = NULL;
}
