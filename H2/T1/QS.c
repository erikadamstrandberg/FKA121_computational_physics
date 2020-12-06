// H2b Variational Monto Carlo
// 
// T1 benchmarking

// Standard C libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// GSL for random number generation
#include <gsl/gsl_rng.h>


// Includes from objective files
#include "init_gsl.h"
void weight(double *w, double *alpha, 
            double *x1, double *y1, double *z1, double *x2, double *y2, double *z2)
{
    double *r1 = malloc(sizeof(double));
    double *r2 = malloc(sizeof(double));
    double *r12 = malloc(sizeof(double));
    
    *r1  = sqrt(pow(*x1,2)       + pow(*y1,2)       + pow(*z1,2));
    *r2  = sqrt(pow(*x2,2)       + pow(*y2,2)       + pow(*z2,2));
    *r12 = sqrt(pow((*x2-*x1),2) + pow((*y2-*y1),2) + pow((*z2-*z1),2));
    
    *w = pow(exp(-2.0*(*r1))*exp(-2.0*(*r2))*exp((*r12)/(2.0*(1.0+*alpha*(*r12)))), 2);

    free(r1);
    free(r2);
    free(r12);
    r1 = NULL;
    r2 = NULL;
    r12 = NULL;

}
// Main 
int main()
{
    // Initializes GSL random number generation
    gsl_rng* gsl_rand = init_gsl();

    // Pick a random number
    double rand;
    rand = gsl_rng_uniform(gsl_rand);

    double alpha = 0.1;
    int N = 10000;
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
    
    FILE *mc = fopen("markov_chain.csv", "w");
    for(int i = 0; i < N; i++)
    {
        fprintf(mc, "%f,%f,%f,%f,%f,%f\n", x1[i], y1[i], z1[i], x2[i], y2[i], z2[i]);
    }
    fclose(mc);
}


