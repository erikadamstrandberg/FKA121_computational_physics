#include <math.h>

void calc_acc(double *a, double *u, double *m, double kappa, double alpha, int size_of_u){
    int i;

    /* Calculating the acceleration on the boundaries */
    a[0] = kappa*(-2*u[0] + u[1])/m[0]
          -alpha*(pow(u[0],2)-pow(u[1]-u[0],2))/m[0];
    a[size_of_u - 1] = kappa*(u[size_of_u - 2] - 2*u[size_of_u - 1])/m[size_of_u - 1]
          -alpha*(pow(u[size_of_u-1]-u[size_of_u-2],2)-pow(u[size_of_u-1],2))/m[size_of_u-1];

    /* Calculating the acceleration of the inner points */
    for (i = 1; i < size_of_u - 1; i++){
        a[i] = kappa*(u[i-1]-2*u[i]+u[i+1])/m[i]
              -alpha*(pow(u[i]-u[i-1],2)-pow(u[i+1]-u[i],2))/m[i];
    }
}
