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
#include "write_to_file.h"

#define NDIM 3


/* Main program */
int main(){
    
    int N = 4;
    int unit_cell = N*N*N;
    int n_atoms = 4*unit_cell;

    double a0_start = cbrt(64);
    double a0_stop  = cbrt(68);
    double lattice_step = 0.005;
    
    int n_a0 = (a0_stop - a0_start)/lattice_step;
    double a0[n_a0];
    double energy[n_a0];
    double L;
    double pos[n_atoms][NDIM];

    for (int i = 0; i < n_a0; i++){
        a0[i] = a0_start + i*lattice_step;
        L = N*a0[i];
        init_fcc(pos, N, a0[i]);
        energy[i] = get_energy_AL(pos, L, n_atoms)/unit_cell;
    }

    FILE *fp = fopen("lattice_pot_energy.csv", "w");
    fprintf(fp, "lattice_const,E_pot\n");

    for(int i = 0; i < n_a0; i++){
	    fprintf(fp, "%f,%f\n", a0[i], energy[i]); 
    }
    fclose(fp);

    print_pos(pos, n_atoms, "initial_AL");
}
