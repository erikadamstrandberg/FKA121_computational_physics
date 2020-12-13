// Standard C libraries
#include <time.h>

// GSL for random number generation
#include <gsl/gsl_rng.h>

// Initializer for GSL random numbers
gsl_rng* init_gsl(){

    const gsl_rng_type *T;
    gsl_rng *Q;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    Q = gsl_rng_alloc(T);
    gsl_rng_set(Q, time(NULL));

    return Q;
}
