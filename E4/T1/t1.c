// Main function to solve E4
//
// Task 1. Implementing the BD3 algorithm
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

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
    double kb = 1.38064852e-23;             // [J/K]
    double rho = 2.65e3;                    // [kg/m^3]
    double w0 = 3.1e3*2.0*PI;             // [Hz]
    double d = 2.79e-6;                      // [m]
    double m = rho*4.0*PI*pow(d/2.0,3)/3.0;
    double temp = 297;                      // [K]
    double vth = sqrt(kb*temp/m);
    // rescaling: time [ms], length [micro m], mass [micro g]
    vth = vth*1e3;
    w0 = w0*1e-3;

    // 
    double mu_low = 1.0/147.3e-3;
    double mu_high = 1.0/48.5e-3;


    double T =  1;                // total simulation time [ms]
    double dt = 0.001;            // timestep [ms] 
    int N = (int) (T/dt);
    int nburn =  (int) (N/10);   
    double c0 = exp(-mu_high*dt);    
    
    double x[N+1];    
    double v[N+1];
    
    x[0] = 0.08;
    v[0] = 0;

    double a = -pow(w0,2)*x[0];
    
    // run burn in 
    for(int t = 1; t < nburn+1; t++){

        v[t] = (dt/2.0)*a+sqrt(c0)*v[t-1]+vth*sqrt(1-c0)*gsl_ran_gaussian(Q,1);
        x[t] = x[t-1] + v[t]*dt;
        a = -pow(w0,2)*x[t];
        v[t] = (dt*sqrt(c0)/2.0)*a+sqrt(c0)*v[t]+vth*sqrt(1-c0)*gsl_ran_gaussian(Q,1);
        
    }

    x[0] = x[nburn];
    v[0] = v[nburn];
    
    // run production
    for(int t = 1; t < N+1; t++){

        v[t] = (dt/2.0)*a+sqrt(c0)*v[t-1]+vth*sqrt(1-c0)*gsl_ran_gaussian(Q,1);
        x[t] = x[t-1] + v[t]*dt;
        a = -pow(w0,2)*x[t];
        v[t] = (dt*sqrt(c0)/2.0)*a+sqrt(c0)*v[t]+vth*sqrt(1-c0)*gsl_ran_gaussian(Q,1);
        
    }
        


    FILE *fp = fopen("timetrail.csv", "w");
    for(int i = 0; i < N+1; i++){
        fprintf(fp, "%f,%f,%f\n", x[i], v[i], dt*i);
    }
    fclose(fp);

}
