// Includes from standard lib
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define PI 3.14159265359


double weight(double *x, double *y, double *z){
    return exp(-pow(*x,2) - pow(*y,2) - pow(*z,2))/pow(PI, (3.0/2.0));
}

double integrand(double *x, double *y, double *z){
    return pow(*x,2) + pow(*x,2)*pow(*y,2) + pow(*x,2)*pow(*y,2)*pow(*z,2);
}

// Main 
int main(){
    int N = 1e8;
    int N_burn = 1e2;

    double x;
    double y;
    double z;
    
    double f;
    double f_trial;
    double x_trial;
    double y_trial;
    double z_trial;
    double r_x;
    double r_y;
    double r_z;
    double delta = 3.1;

    double accept_rand;

    srand(time(NULL));
    x = (double) rand()/(double) RAND_MAX;
    y = (double) rand()/(double) RAND_MAX;
    z = (double) rand()/(double) RAND_MAX;

    double integral = 0.0;

    int accept = 0;
    for(int i = 0; i < N-1; i++){
        r_x = (double) rand()/(double) RAND_MAX;
        r_y = (double) rand()/(double) RAND_MAX;
        r_z = (double) rand()/(double) RAND_MAX;

        x_trial = x + delta*(r_x-0.5);
        y_trial = y + delta*(r_y-0.5);
        z_trial = z + delta*(r_z-0.5);

        f = weight(&x, &y, &z);
        f_trial = weight(&x_trial, &y_trial, &z_trial);

        accept_rand = (double) rand()/(double) RAND_MAX;
        
        if(f_trial > f || f_trial/f > accept_rand){
            x = x_trial;    
            y = y_trial;    
            z = z_trial;
            accept += 1;
        }

        if(i > N_burn){
            integral += integrand(&x, &y, &z); 
        }
    }

    integral += integrand(&x, &y, &z);
    integral = integral/(N-N_burn);

    printf("Exact integral: %f\n", 7.0/8.0);
    printf("Calculated integral: %f\n", integral);
    printf("Acceptance ratio: %f\n\n", (double) accept/(double) N);
}
