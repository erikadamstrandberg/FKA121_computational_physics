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

void mean(double *mean, double *array, int start, int stopp);
void sigma2(double *sigma2, double *array, int start, int stopp);
void correlation_function(double *phi, double *E_l, double *sigma2_E_l, int N, int k);
void block_averaging(double *ns, double *E_l, double *sigma2_E_l, int N, int B);

// Main 
int main()
{
    // Initializes GSL random number generation
    gsl_rng* gsl_rand = init_gsl();

    // Pick a random number
    double rand;
    rand = gsl_rng_uniform(gsl_rand);

    double alpha = 0.1;
    int N_tot = 50000;
    int burn_in = 2000;
    int N = N_tot - burn_in;

    double x1[N_tot];
    double y1[N_tot];
    double z1[N_tot];
    double x1_t;
    double y1_t;
    double z1_t;
    
    double x2[N_tot];
    double y2[N_tot];
    double z2[N_tot];
    double x2_t;
    double y2_t;
    double z2_t;

    double w;
    double w_t;
    double delta = 1.0;
    int accept = 0;

    double E_l[N];

    x1[0] = delta*(gsl_rng_uniform(gsl_rand)-0.5);
    y1[0] = delta*(gsl_rng_uniform(gsl_rand)-0.5);
    z1[0] = delta*(gsl_rng_uniform(gsl_rand)-0.5);

    x2[0] = delta*(gsl_rng_uniform(gsl_rand)-0.5);
    y2[0] = delta*(gsl_rng_uniform(gsl_rand)-0.5);
    z2[0] = delta*(gsl_rng_uniform(gsl_rand)-0.5);

    for(int i = 0; i < N_tot; i++)
    {
        x1_t = x1[i] + delta*(gsl_rng_uniform(gsl_rand)-0.5);
        y1_t = y1[i] + delta*(gsl_rng_uniform(gsl_rand)-0.5);
        z1_t = z1[i] + delta*(gsl_rng_uniform(gsl_rand)-0.5);
        x2_t = x2[i] + delta*(gsl_rng_uniform(gsl_rand)-0.5);
        y2_t = y2[i] + delta*(gsl_rng_uniform(gsl_rand)-0.5);
        z2_t = z2[i] + delta*(gsl_rng_uniform(gsl_rand)-0.5);

        weight(&w, &alpha, &x1[i], &y1[i], &z1[i], &x2[i], &y2[i], &z2[i]); 
        weight(&w_t, &alpha, &x1_t, &y1_t, &z1_t, &x2_t, &y2_t, &z2_t);

        if(w_t > w || w_t/w > gsl_rng_uniform(gsl_rand))
        {
           x1[i+1] = x1_t; 
           y1[i+1] = y1_t; 
           z1[i+1] = z1_t; 
           x2[i+1] = x2_t; 
           y2[i+1] = y2_t; 
           z2[i+1] = z2_t;
           accept += 1;
        } 
        else 
        {
           x1[i+1] = x1[i]; 
           y1[i+1] = y1[i]; 
           z1[i+1] = z1[i]; 
           x2[i+1] = x2[i]; 
           y2[i+1] = y2[i]; 
           z2[i+1] = z2[i];
        }
    }

    printf("Acceptance ratio: %f\n", (double) accept/((double) N));
    

    local_energy(E_l, &alpha, N_tot, burn_in, x1, y1, z1, x2, y2, z2);
    print_markov(N_tot, x1, y1, z1, x2, y2, z2, "markov_chain");
    print_1d_array(E_l, N, "local_energy");
    
    double mean_E_l = 0.0;
    double sigma2_E_l = 0.0;
    mean(&mean_E_l, E_l, 0, N);
    sigma2(&sigma2_E_l, E_l, 0, N);

    int k = 30;
    double phi[k];    
    correlation_function(phi, E_l, &sigma2_E_l, N, k);
    print_1d_array(phi, k, "correlation");

    int B = 120;
    double ns_b[B];
    block_averaging(ns_b, E_l, &sigma2_E_l, N, B);
    print_1d_array(ns_b, B, "block_average");

    double ns;
    double limit = 0.135;
    int p = 0;
    while (phi[p] > limit)
    {
        p += 1;
    }

    printf("Mean of E_l: %f\n", mean_E_l);
    printf("Variance of E_l: %f\n", sigma2_E_l);
    printf("Statistical inefficiency: %d\n\n", p);
    printf("Estimate: %f +- %f\n",  mean_E_l, sqrt(sigma2_E_l)/(sqrt((double) N/(double) p)));
}

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
    
    free(mean_array);
    mean_array = NULL;
}

void correlation_function(double *phi, double *E_l, double *sigma2_E_l, int N, int k)
{
    double mean_forward;
    double mean_lag;

    for(int i = 0; i < k; i++)
    {
        mean(&mean_forward, E_l, i, N);
        mean(&mean_lag, E_l, 0, N-i);
        for(int j = 0; j < N-i; j++)
        {
            phi[i] += (E_l[i+j] - mean_forward)*(E_l[j] - mean_lag);
        }
        phi[i] = phi[i]/(*sigma2_E_l*(N-i));
    }
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
