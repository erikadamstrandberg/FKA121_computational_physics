#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// GSL for random number generation
#include <gsl/gsl_rng.h>

void local_energy(double *E_l, double *alpha, 
                  double *x1, double *y1, double *z1, double *x2, double *y2, double *z2)
{
    // Allocating variables
    double *r1 = malloc(sizeof(double));
    double *r2 = malloc(sizeof(double));
    double *r12 = malloc(sizeof(double));
    double *snd_term = malloc(sizeof(double));

    // Spliting up the calculation into separate pieces
    *r1  = sqrt(pow(*x1,2)         + pow(*y1,2)         + pow(*z1,2));
    *r2  = sqrt(pow(*x2,2)         + pow(*y2,2)         + pow(*z2,2));
    *r12 = sqrt(pow((*x2-*x1),2) + pow((*y2-*y1),2) + pow((*z2-*z1),2));

    *snd_term = (((*r2)*(*x1)-(*r1)*(*x2))*((*x1)-(*x2))+
                 ((*r2)*(*y1)-(*r1)*(*y2))*((*y1)-(*y2))+
                 ((*r2)*(*z1)-(*r1)*(*z2))*((*z1)-(*z2)))
                 /((*r1)*(*r2)*(*r12)*pow((1+(*alpha)*(*r12)),2));

    // Actual calculation of the local energy   &E_l??
    *E_l = -4.0 + (*snd_term) -1.0/((*r12)*pow((1.0 + (*alpha)*(*r12)),3)) 
           -1.0/(4.0*pow((1.0 + (*alpha)*(*r12)),4)) + 1.0/(*r12); 


    // Free all variables
    free(r1); r1 = NULL;
    free(r2); r2 = NULL;
    free(r12); r12 = NULL;
    free(snd_term); snd_term = NULL;

}

void weight(double *w, double *alpha, 
            double *x1, double *y1, double *z1, double *x2, double *y2, double *z2)
{
    // Allocating variables
    double *r1 = malloc(sizeof(double));
    double *r2 = malloc(sizeof(double));
    double *r12 = malloc(sizeof(double));
    
    // Splitting the wavefunction into readable pieces
    *r1  = sqrt(pow(*x1,2)       + pow(*y1,2)       + pow(*z1,2));
    *r2  = sqrt(pow(*x2,2)       + pow(*y2,2)       + pow(*z2,2));
    *r12 = sqrt(pow((*x2-*x1),2) + pow((*y2-*y1),2) + pow((*z2-*z1),2));
    
    // Actual calculation of the weight function
    *w = pow(exp(-2.0*(*r1))*exp(-2.0*(*r2))*exp((*r12)/(2.0*(1.0+*alpha*(*r12)))), 2);

    // Free all variables
    free(r1); r1 = NULL;
    free(r2); r2 = NULL;
    free(r12); r12 = NULL;
    
}

void metropolis_move(double *x1, double *y1, double *z1, 
                     double *x2, double *y2, double *z2,
                     double *delta,
                     double *alpha,
                     int *accept,
                     gsl_rng *gsl_rand)
{
    // Allocating variables
    double *x1_t = malloc(sizeof(double));
    double *y1_t = malloc(sizeof(double));
    double *z1_t = malloc(sizeof(double));
    double *x2_t = malloc(sizeof(double));
    double *y2_t = malloc(sizeof(double));
    double *z2_t = malloc(sizeof(double));

    // Variables for storing the weight function
    double *w = malloc(sizeof(double));
    double *w_t = malloc(sizeof(double));

    // Generating the trial step
    *x1_t = *x1 + (*delta)*(gsl_rng_uniform(gsl_rand)-0.5);
    *y1_t = *y1 + (*delta)*(gsl_rng_uniform(gsl_rand)-0.5);
    *z1_t = *z1 + (*delta)*(gsl_rng_uniform(gsl_rand)-0.5);
    *x2_t = *x2 + (*delta)*(gsl_rng_uniform(gsl_rand)-0.5);
    *y2_t = *y2 + (*delta)*(gsl_rng_uniform(gsl_rand)-0.5);
    *z2_t = *z2 + (*delta)*(gsl_rng_uniform(gsl_rand)-0.5);

    // Calculating the proability distribution for last step
    // and the trial step
    weight(w, alpha, x1, y1, z1, x2, y2, z2);  
    weight(w_t, alpha, x1_t, y1_t, z1_t, x2_t, y2_t, z2_t);

    // Metropolis algorithm
    if((*w_t) > (*w) || (*w_t)/(*w) > gsl_rng_uniform(gsl_rand))
    {
       // Accepting the trial step
       *x1 = *x1_t; 
       *y1 = *y1_t; 
       *z1 = *z1_t; 
       *x2 = *x2_t; 
       *y2 = *y2_t; 
       *z2 = *z2_t;
       *accept += 1;
    }
    
    // Freeing all variables
    free(x1_t); x1_t = NULL; 
    free(y1_t); y1_t = NULL; 
    free(z1_t); z1_t = NULL; 
    free(x2_t); x2_t = NULL; 
    free(y2_t); y2_t = NULL; 
    free(z2_t); z2_t = NULL;
    free(w); w = NULL;
    free(w_t); w_t = NULL;    

}

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
