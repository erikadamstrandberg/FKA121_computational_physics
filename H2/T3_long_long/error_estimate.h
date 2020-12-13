void mean(double *mean, double *array, int start, int stopp);
void sigma2(double *sigma2, double *array, int start, int stopp);
void correlation_function(double *phi, double *E_l, double *sigma2_E_l, int N, int k);
void block_averaging(double *ns, double *E_l, double *sigma2_E_l, int N, int B);
void estimate_ns(double *ns, double *E_l, double sigma2_E_l, int N);
void initialize_mcmc(double *ns, double *alpha,
            double *E_l_mean, double *E_l2_mean,
            double *fst_term, double *snd_term,
            double *delta, int *accept, int N_ns, int burn_in, 
            double *x1, double *y1, double *z1, double *x2, double *y2, double *z2,
            gsl_rng *gsl_rand);

