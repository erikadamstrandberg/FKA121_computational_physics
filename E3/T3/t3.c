// Includes from standard lib
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define PI 3.14159265359


double weight(double *x, double *y, double *z){
    return pow(PI, -3.0/2.0)*exp(-(pow(*x,2)+pow(*y,2)+pow(*z,2)));
}

double integrand(double *x, double *y, double *z){
    return pow(*x,2) + pow(*x,2)*pow(*y,2) + pow(*x,2)*pow(*y,2)*pow(*z,2);
}

// Main 
int main(){
    int N = 100000;
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
    double delta = 2.0;

    srand(time(NULL));
    x[0] = (double) rand()/(double) RAND_MAX;
    y[0] = (double) rand()/(double) RAND_MAX;
    z[0] = (double) rand()/(double) RAND_MAX;

    int accept = 0;
    for(int i = 0; i < N-1; i++){
        r_x = (double) rand()/(double) RAND_MAX;
        r_y = (double) rand()/(double) RAND_MAX;
        r_z = (double) rand()/(double) RAND_MAX;
        x_trial = x[i] + delta*(r_x-0.5);
        y_trial = y[i] + delta*(r_y-0.5);
        z_trial = z[i] + delta*(r_z-0.5);

        f = weight(&x[i], &y[i], &z[i]);
        f_trial = weight(&x_trial, &y_trial, &z_trial);

        if(f_trial > f || f_trial/f > r_x){
            x[i+1] = x_trial;    
            y[i+1] = y_trial;    
            z[i+1] = z_trial;

            accept += 1;   

        }else{
            x[i+1] = x[i];    
            y[i+1] = y[i];    
            z[i+1] = z[i];    
        }
    }


    double integral = 0.0;
    for(int i = N_burn; i < N; i++){
        integral += integrand(&x[i], &y[i], &z[i]); 
    }


    integral = integral/(N-N_burn);
    printf("%f\n", integral);
    printf("%f\n", (double) accept/(double) N);
}
