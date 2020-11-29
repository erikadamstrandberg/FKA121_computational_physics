#define NDIM 3

#ifndef _alpotential_h
#define _alpotential_h

extern void get_forces_AL(double[][3] , double[][3], double, int);
extern double get_energy_AL(double[][3], double, int);
extern double get_virial_AL(double[][3], double, int);
extern void random_displacement(double pos[][NDIM], int n_atoms, double interval);
extern void verlet_timestep(double pos[][NDIM], double v[][NDIM], double f[][NDIM], 
                     int n_atoms, 
                     double dt, 
                     double m, 
                     double L);

extern double get_kinetic_energy(double v[][NDIM], int n_atoms, double m);

#endif
