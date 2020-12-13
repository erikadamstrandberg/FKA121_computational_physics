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

void estimate_ns(double *ns, double *E_l, double sigma2_E_l, int N)
{

    // Caculate and writes the correlation function
    int k = 30;
    double phi[k];    
    correlation_function(phi, E_l, &sigma2_E_l, N, k);

    // Calculates and writes the block averaging
    int B = 800;
    double ns_b[B];
    block_averaging(ns_b, E_l, &sigma2_E_l, N, B);

    // Estimate ns from correlation function
    double limit = 0.135;
    double ns_from_corr;
    int index_after = 0;
    while (phi[index_after] > limit)
    {
        index_after += 1;
    }

    ns_from_corr = index_after*(phi[index_after-1]-limit)/(phi[index_after-1]-phi[index_after]) + 
                  (index_after-1)*(limit-phi[index_after])/(phi[index_after-1]-phi[index_after]); 

    // Estimate from block average
    double ns_from_block = 0;
    int blocks_from_end = 200;
    for(int i = B-blocks_from_end; i < B; i++)
    {
        ns_from_block += ns_b[i];
    }
    ns_from_block = ns_from_block/blocks_from_end;

    if (ns_from_block > ns_from_corr)
    {
        *ns = ns_from_block;
    }
    else 
    {
        *ns = ns_from_corr;
    }
}
