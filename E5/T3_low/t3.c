// Main function to solve E4
//
// Task 1. Implementing the BD3 algorithm
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <gsl/gsl_randist.h>

#define PI 3.141592653589

// Main 
int main(){
    
    // setup  for random number generation    
    const gsl_rng_type *U;
    gsl_rng *Q;
    gsl_rng_env_setup();
    U = gsl_rng_default;
    Q = gsl_rng_alloc(U);
    gsl_rng_set(Q, time(NULL));
    
    // define constants
    double kb = 1.38064852e-23;              // [J/K]
    double rho = 2.65e3;                     // [kg/m^3]
    double w0 = 3.1e3*2.0*PI;                // [Hz]
    double d = 2.79e-6;                      // [m]
    double m = rho*4.0*PI*pow(d/2.0,3)/3.0;
    double temp = 297;                       // [K]
    double vth = sqrt(kb*temp/m);
    // rescaling: time [ms], length [micro m], mass [micro g]
    vth = vth*1e3;
    w0 = w0*1e-3;

    // 
    double mu_low = 1.0/147.3e-3;
    double mu_high = 1.0/48.5e-3;


    double T = 5.0;                // total simulation time [ms]
    double dt = 0.001;            // timestep [ms] 
    int N = (int) (T/dt);
    double c0 = exp(-mu_low*dt);
    
    double x[N+1];    
    double v[N+1];
   
    char filename[] = "data/timetrail";
    char run_buffer[10];

    int how_many_runs = 10000;
    for(int p = 0; p < how_many_runs; p++){
        // run production
        

        x[0] = 0.1;
        v[0] = 2.0;

        double a = -pow(w0,2)*x[0];
        
        for(int t = 1; t < N+1; t++){
            v[t] = (dt/2.0)*a+sqrt(c0)*v[t-1]+vth*sqrt(1-c0)*gsl_ran_gaussian(Q,1);
            x[t] = x[t-1] + v[t]*dt;
            a = -pow(w0,2)*x[t];
            v[t] = (dt*sqrt(c0)/2.0)*a+sqrt(c0)*v[t]+vth*sqrt(1-c0)*gsl_ran_gaussian(Q,1);
        
        }


        char *filename_csv = malloc(40*sizeof(char));
        sprintf(run_buffer, "%d", p);
        strcpy(filename_csv, filename);
        strcat(filename_csv, run_buffer);
        strcat(filename_csv, ".csv");
        
        FILE *fp = fopen(filename_csv, "w");
        for(int i = 0; i < N+1; i++){
            fprintf(fp, "%f,%f,%f\n", x[i], v[i], dt*i);
        }
        fclose(fp);
        free(filename_csv); filename_csv = NULL;
    }
}
