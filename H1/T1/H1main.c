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
    
    /*
     Code for generating a uniform random number between 0 and 1. srand should only
     be called once.
    */
    /*
     srand(time(NULL));
     double random_value;
     random_value = (double) rand() / (double) RAND_MAX;
    */
    
    /*
     Descriptions of the different functions in the files H1lattice.c and
     H1potential.c are listed below.
    */
    
    /* 
     Function that generates a fcc lattice in units of [Å]. Nc is the number of 
     primitive cells in each direction and a0 is the lattice parameter. The
     positions of all the atoms are stored in pos which should be a matrix of the
     size N x 3, where N is the number of atoms. The first, second and third column
     correspond to the x,y and z coordinate respectively.
    */
    double a0 = 4.05;  
    int N = 4;
    int n_atoms = 4*N*N*N;

    double pos[n_atoms][3];
    init_fcc(pos, N, a0);
    
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
