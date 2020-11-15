/*
 H1main.c
 
 Created by Anders Lindman on 2013-10-31.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "H1lattice.h"
#include "H1potential.h"


/* Main program */
int main()
{
    int ndim = 3;
    double a0 = 4.05;  
    int N = 4;
    int n_atoms = 4*N*N*N;
    double pos[n_atoms][ndim];
    double rand_interval = 0.065*a0;

    init_fcc(pos, N, a0);

    srand(time(NULL));
    double random_value;
    double random_disp;

    for(int i = 0; i < n_atoms; i++){
        for(int j = 0; j < ndim; j++){
            random_value = (double) rand() / (double) RAND_MAX;
            random_disp = random_value*2*rand_interval-rand_interval;
            pos[i][j] += random_disp;
            printf("%f\n", pos[i][j]);
        }
    }
    
    
    /*
    for(int i = 0; i < n_atoms; i++){
        printf("%f\n", pos[i][0]);
    }
    */

    /* 
     Function that calculates the potential energy in units of [eV]. pos should be
     a matrix containing the positions of all the atoms, L is the length of the 
     supercell and N is the number of atoms.
    */
     double L = N*a0;
     double energy;
     energy = get_energy_AL(pos, L, n_atoms);
    
     printf("%f\n", energy);
    /* 
     Function that calculates the virial in units of [eV]. pos should be a matrix
     containing the positions of all the atoms, L is the length of the supercell 
     and N is the number of atoms.
    */
    /*
     double virial;
     virial = get_virial_AL(pos, L, N);
    */
    
    /*
     Function that calculates the forces on all atoms in units of [eV/Å]. the 
     forces are stored in f which should be a matrix of size N x 3, where N is the
     number of atoms and column 1,2 and 3 correspond to the x,y and z component of
     the force resepctively . pos should be a matrix containing the positions of 
     all the atoms, L is the length of the supercell and N is the number of atoms.
    */
    /*
     get_forces_AL(f,pos, L, N);
    */
    
    
    
}
