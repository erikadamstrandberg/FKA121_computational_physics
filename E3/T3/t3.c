// Includes from standard lib
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <gsl/gsl_rng.h>

#define PI 3.14159265359


double weight(double *x, double *y, double *z){
    return exp(-pow(*x,2) - pow(*y,2) - pow(*z,2))/pow(PI, (3.0/2.0));
}

double integrand(double *x, double *y, double *z){
    return pow(*x,2) + pow(*x,2)*pow(*y,2) + pow(*x,2)*pow(*y,2)*pow(*z,2);
}

// Main 
int main(){
    int N = 3e4;
    int N_burn = 1000;

    double x[N];
    double y[N];
    double z[N];
    
    double f;
    double f_trial;
    double x_trial;
    double y_trial;
    double z_trial;
    double r_x;
    double r_y;
    double r_z;
    double delta = 2.5;

    double accept_rand;
    
    
    const gsl_rng_type *T;
    gsl_rng *Q;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    Q = gsl_rng_alloc(T);
    gsl_rng_set(Q, time(NULL));

    x[0] = gsl_rng_uniform(Q);
    y[0] = gsl_rng_uniform(Q);
    z[0] = gsl_rng_uniform(Q);

    int accept = 0;
    for(int i = 0; i < N-1; i++){
        r_x = gsl_rng_uniform(Q);
        r_y = gsl_rng_uniform(Q);
        r_z = gsl_rng_uniform(Q);
        
        x_trial = x[i] + delta*(r_x-0.5);
        y_trial = y[i] + delta*(r_y-0.5);
        z_trial = z[i] + delta*(r_z-0.5);

        f = weight(&x[i], &y[i], &z[i]);
        f_trial = weight(&x_trial, &y_trial, &z_trial);

        accept_rand = gsl_rng_uniform(Q);
        
        if(f_trial > f || f_trial/f > accept_rand){
            x[i+1] = x_trial;    
            y[i+1] = y_trial;    
            z[i+1] = z_trial;
            accept += 1;
        } else { 
            x[i+1] = x[i];    
            y[i+1] = y[i];    
            z[i+1] = z[i];
        }
    }


    FILE *fp = fopen("walkers.csv", "w");
    fprintf(fp, "x,y,z,weight\n");
    for(int i = 0; i < N; i++){
        fprintf(fp, "%f,%f,%f,%f\n", x[i], y[i], z[i], weight(&x[i], &y[i], &z[i]));
    }
    fclose(fp);

    FILE *fn = fopen("N.csv", "w");
    fprintf(fp, "N\n");
    fprintf(fp, "%d,%d", N, N);
    fclose(fn);

    double integral = 0.0;
    for(int i = N_burn; i < N; i++){
        integral += integrand(&x[i], &y[i], &z[i]); 
    }

    integral = integral/(N-N_burn);

    printf("Exact integral: %f\n", 7.0/8.0);
    printf("Calculated integral: %f\n", integral);
    printf("Acceptance ratio: %f\n\n", (double) accept/(double) N);
}
