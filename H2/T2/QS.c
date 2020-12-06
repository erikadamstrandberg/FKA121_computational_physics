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

// Main 
int main()
{
    // Initializes GSL random number generation
    gsl_rng* gsl_rand = init_gsl();

    // Pick a random number
    double rand;
    rand = gsl_rng_uniform(gsl_rand);

    double alpha = 0.1;
    int N = 20000;
    double x1[N];
    double y1[N];
    double z1[N];
    double x1_t;
    double y1_t;
    double z1_t;
    
    double x2[N];
    double y2[N];
    double z2[N];
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


    for(int i = 0; i < N; i++)
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

    printf("%f\n", (double) accept/((double) N));
    
    local_energy(E_l, &alpha, N, x1, y1, z1, x2, y2, z2);
    print_markov(E_l, N, x1, y1, z1, x2, y2, z2, "markov_chaimarkov_chain");
}
