// Includes from standard lib
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define PI 3.14159265359


double weight(double *x, double *y, double *z){
    return pow(PI, -3.0/2.0)*exp(-(pow(*x,2) + pow(*y,2) + pow(*z,2)));
}

double integrand(double *x, double *y, double *z){
    return pow(*x,2) + pow(*x,2)*pow(*y,2) + pow(*x,2)*pow(*y,2)*pow(*z,2);
}

// Main 
int main(){
    int N = 300000;
    int N_burn = 20000;

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
    double delta = 2.4;

    double accept_rand;

    srand(time(NULL));
    x[0] = (double) rand()/(double) RAND_MAX;
    y[0] = (double) rand()/(double) RAND_MAX;
    z[0] = (double) rand()/(double) RAND_MAX;

    int accept = 0;
    int non_accept = 0;
    for(int i = 0; i < N-1; i++){
        r_x = (double) rand()/(double) RAND_MAX;
        r_y = (double) rand()/(double) RAND_MAX;
        r_z = (double) rand()/(double) RAND_MAX;

        x_trial = x[i] + delta*(r_x-0.5);
        y_trial = y[i] + delta*(r_y-0.5);
        z_trial = z[i] + delta*(r_z-0.5);

        f = weight(&x[i], &y[i], &z[i]);
        f_trial = weight(&x_trial, &y_trial, &z_trial);

        accept_rand = (double) rand()/(double) RAND_MAX;
        
        if(f_trial > f){
            x[i+1] = x_trial;    
            y[i+1] = y_trial;    
            z[i+1] = z_trial;
        } else {
            if (f_trial/f > accept_rand){
                x[i+1] = x_trial;    
                y[i+1] = y_trial;    
                z[i+1] = z_trial;
                accept += 1;
            } else {
                x[i+1] = x[i];    
                y[i+1] = y[i];    
                z[i+1] = z[i];
                non_accept += 1;
            }
        }
    }


    FILE *fp = fopen("walkers.csv", "w");
    fprintf(fp, "x,y,z\n");
    for(int i = 0; i < N; i++){
        fprintf(fp, "%f,%f,%f\n", x[i], y[i], z[i]);
    }
    fclose(fp);

    double integral = 0.0;
    for(int i = N_burn; i < N; i++){
        integral += integrand(&x[i], &y[i], &z[i]); 
    }


    integral = integral/(N-N_burn);

    printf("Exact integral: %f\n", 7.0/8.0);
    printf("Calculated integral: %f\n", integral);
    printf("Acceptance ratio: %f\n\n", (double) accept/(double) non_accept);
}
