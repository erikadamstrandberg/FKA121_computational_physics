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

#%%

U_tot_mean = sum(U_tot)/len(U_tot)
print(U_tot_mean)

#%%

mode_energy = np.genfromtxt('energy_norm.csv', delimiter=',', skip_header=1)
timesteps = len(mode_energy[0,:])

fig, ax = plt.subplots()

ax.plot(time, mode_energy[:, 0], color='black', linewidth=3, label=r'$E_1(t)$')
ax.plot(time, mode_energy[:, 1], label=r'$E_2(t)$',linewidth=10)
ax.plot(time, mode_energy[:, 2], label=r'$E_3(t)$',linewidth=5)
ax.plot(time, mode_energy[:, 3], label=r'$E_4(t)$',linewidth=2)
ax.plot(time, mode_energy[:, 4], label=r'$E_5(t)$',linewidth=0.5)

ax.set_xlabel('Energy (arb.unit)')
ax.set_ylabel('Time (arb.unit)')  
ax.grid()
ax.legend(fontsize=16)

plt.savefig('mode_ortogonality.pdf', format='pdf', bbox_inches='tight')

#%%


