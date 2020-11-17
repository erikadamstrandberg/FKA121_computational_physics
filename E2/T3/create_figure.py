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
time  = np.log10(array[:, 2])

fig, ax = plt.subplots()
ax.plot(time, U_kin)
ax.plot(time, U_pot)
ax.plot(time, U_tot)

#ax.set_title(r'Conservation of energy, $\alpha=0.01$', fontsize='20')
ax.set_xlabel(r'Time', fontsize='16')
ax.set_ylabel(r'Energy', fontsize='16')
ax.grid()
plt.savefig("energy.pdf", format="pdf", bbox_inches="tight")

U_tot_mean = sum(U_tot)/len(U_tot)
print(U_tot_mean)


mode_energy = np.genfromtxt('energy_norm.csv', delimiter=',', skip_header=1)
n_timesteps = len(mode_energy[:,0])
n_modes = len(mode_energy[0,:])

E_average = np.zeros([n_timesteps, n_modes])
for t in range(n_timesteps):
    for i in range(n_modes):
        E_average[t][i] = np.sum(mode_energy[0:t,i])/(t+1)
    
E_average = np.log10(E_average)
fig, ax = plt.subplots()
for i in range(n_modes):
    ax.plot(time, E_average[:,i])
ax.set_ylabel('Average Energy (arb.unit)')
ax.set_xlabel('time (arb.units)')  
ax.grid()
ax.legend(fontsize=16)
plt.savefig("average_energy.pdf", format="pdf", bbox_inches="tight")

fig, ax = plt.subplots()
mode_energy = np.log10(mode_energy)
ax.plot(time, mode_energy[:, 0])
ax.plot(time, mode_energy[:, 1])
ax.plot(time, mode_energy[:, 2])
ax.plot(time, mode_energy[:, 3])
ax.plot(time, mode_energy[:, 5])


ax.set_xlabel('Time (arb.unit)')  
ax.set_ylabel('Energy (arb.unit)')
ax.grid()
ax.legend(fontsize=16)
plt.savefig("energy_modes.pdf", format="pdf", bbox_inches="tight")
plt.show()

