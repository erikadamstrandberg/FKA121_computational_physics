// H3b Time independent quantum mechanics 
// 
// T1 making the FFT work!

// Standard C libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include <fftw3.h>

// Our own includes
#include "fft.h"
#include "write_to_file.h"

// Constants
#define PI 3.14159

// Main 
int main()
{
    int N = 1000;
    double x[N];
    complex x_complex[N];
    double t[N];

    double f  = 1;
    double dt = 0.01;
    double a  = 1;

    for(int i = 0; i < N; i++)
    {
        x[i] = a*cos(2*PI*f*i*dt);
        x_complex[i] = a*cos(2*PI*f*i*dt);
        t[i] = dt*i;
    }

    print_1d_array(x, N, "signal");
    print_1d_array(t, N, "time");

    double spectrum[N];
    double freq[N];

    fftw_complex out[N];
    fftw_plan p;
    p = fftw_plan_dft_1d(N, x_complex, out, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(p);

    complex out_complex[N];
    for(int i = 0; i < N; i++)
    {   
        out_complex[i] = creal(out[i]) + I*cimag(out[i]);
    }
    
    print_complex_array(out_complex, N, "out");
}
