void local_energy(double *E_l, double *alpha, 
                  double *x1, double *y1, double *z1, double *x2, double *y2, double *z2);
void weight(double *w, double *alpha, 
            double *x1, double *y1, double *z1, double *x2, double *y2, double *z2);
void metropolis_move(double *x1, double *y1, double *z1, 
                     double *x2, double *y2, double *z2,
                     double *delta,
                     double *alpha,
                     int *accept,
                     gsl_rng *gsl_rand);
void grad_alpha_ln_phi(double *grad_ln_phi, double *alpha, 
                       double *x1, double *y1, double *z1, double *x2, double *y2, double *z2);
