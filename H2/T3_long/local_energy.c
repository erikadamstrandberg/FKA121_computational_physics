#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void local_energy(double *E_l, double *alpha, int N, int burn_in, 
                  double *x1, double *y1, double *z1, double *x2, double *y2, double *z2)
{
    // Allocating variables
    double *r1 = malloc(sizeof(double));
    double *r2 = malloc(sizeof(double));
    double *r12 = malloc(sizeof(double));
    double *snd_term = malloc(sizeof(double));

    // Loops over the entire Markov chain
    for (int i = burn_in; i < N; i++)
    {
        // Spliting up the calculation into separate pieces
        *r1  = sqrt(pow(x1[i],2)         + pow(y1[i],2)         + pow(z1[i],2));
        *r2  = sqrt(pow(x2[i],2)         + pow(y2[i],2)         + pow(z2[i],2));
        *r12 = sqrt(pow((x2[i]-x1[i]),2) + pow((y2[i]-y1[i]),2) + pow((z2[i]-z1[i]),2));

        *snd_term = (((*r2)*x1[i]-(*r1)*x2[i])*(x1[i]-x2[i])+
                     ((*r2)*y1[i]-(*r1)*y2[i])*(y1[i]-y2[i])+
                     ((*r2)*z1[i]-(*r1)*z2[i])*(z1[i]-z2[i]))
                    /((*r1)*(*r2)*(*r12)*pow((1+(*alpha)*(*r12)),2));

        // Actual calculation of the local energy
        E_l[i-burn_in] = -4.0 + (*snd_term) -1.0/((*r12)*pow((1.0 + (*alpha)*(*r12)),3)) 
                         -1.0/(4.0*pow((1.0 + (*alpha)*(*r12)),4)) + 1.0/(*r12); 

    }

    // Free all variables
    free(r1);
    free(r2);
    free(r12);
    free(snd_term);
    r1 = NULL;
    r2 = NULL;
    r12 = NULL;
    snd_term = NULL;

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
    free(r1);
    free(r2);
    free(r12);
    r1 = NULL;
    r2 = NULL;
    r12 = NULL;
}
