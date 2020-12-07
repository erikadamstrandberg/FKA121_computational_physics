#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void mean(double *mean, double *array, int start, int stop)
{
    *mean = 0.0;
    for(int i = start; i < stop; i++)
    {
        *mean += array[i];
    }
    *mean = *mean/(stop-start);
}

void sigma2(double *sigma2, double *array, int start, int stop)
{
    double *mean_array = malloc(sizeof(double));
    mean(mean_array, array, start, stop);

    *sigma2 = 0.0;
    for(int i = start; i < stop; i++)
    {
        *sigma2 += pow(array[i], 2);
    }
    *sigma2 = *sigma2/(stop-start) - pow(*mean_array,2);
    
    free(mean_array);
    mean_array = NULL;
}

void correlation_function(double *phi, double *E_l, double *sigma2_E_l, int N, int k)
{
    double mean_forward;
    double mean_lag;

    for(int i = 0; i < k; i++)
    {
        mean(&mean_forward, E_l, i, N);
        mean(&mean_lag, E_l, 0, N-i);
        for(int j = 0; j < N-i; j++)
        {
            phi[i] += (E_l[i+j] - mean_forward)*(E_l[j] - mean_lag);
        }
        phi[i] = phi[i]/(*sigma2_E_l*(N-i));
    }
}

void block_averaging(double *ns, double *E_l, double *sigma2_E_l, int N, int B)
{ 
    double *Mb = malloc(sizeof(double));
    double *sigma2_each_block = malloc(sizeof(double));
    for(int i = 1; i <= B; i++)
    {
        *Mb = N/i;
        double *F = malloc((*Mb)*sizeof(double));
        for(int j = 0; j < *Mb; j++) 
        {
            F[j] = 0.0;
            for(int k = 0; k < i; k++)
            {
                F[j] += E_l[k+j*i]; 
            }
            F[j] = F[j]/i;
        }
        sigma2(sigma2_each_block, F, 0, *Mb);

        ns[i-1] = i*(*sigma2_each_block)/(*sigma2_E_l);
        free(F);
        F = NULL;
    }
    free(Mb);
    free(sigma2_each_block);
    Mb = NULL;
    sigma2_each_block = NULL;
}
