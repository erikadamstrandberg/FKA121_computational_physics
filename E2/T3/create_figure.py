#!/usr/bin/env python
###############################################################################
# Create figure
# plt.savefig('name_of_figure.pdf', bbox_inches='tight')
###############################################################################

import numpy as np
import matplotlib.pyplot as plt


#%%

array = np.genfromtxt('energy_ord.csv', delimiter=',', skip_header=1)
U_kin = array[:, 0]
U_pot = array[:, 1]
U_tot = array[:, 0] + array[:, 1]
time  = array[:, 2]

fig, ax = plt.subplots()
ax.plot(time, U_kin)
ax.plot(time, U_pot)
ax.plot(time, U_tot)

ax.set_title(r'Conservation of energy, $\alpha=0.01$', fontsize='20')
ax.set_xlabel(r'Time', fontsize='16')
ax.set_ylabel(r'Energy', fontsize='16')
ax.grid()

U_tot_mean = sum(U_tot)/len(U_tot)
print(U_tot_mean)


mode_energy = np.genfromtxt('energy_norm.csv', delimiter=',', skip_header=1)
timesteps = len(mode_energy[0,:])

fig, ax = plt.subplots()

ax.plot(np.log10(time), np.log10(mode_energy[:, 0]))
ax.plot(np.log10(time), np.log10(mode_energy[:, 1]))
ax.plot(np.log10(time), np.log10(mode_energy[:, 2]))


ax.set_xlabel('Energy (arb.unit)')
ax.set_ylabel('Time (arb.unit)')  
ax.grid()
ax.legend(fontsize=16)

