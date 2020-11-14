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

ax.set_xlabel('time (arb.unit)')
ax.set_ylabel('signal (arb.unit)')
ax.grid()

ax.set_xlim([0, 200])

#%%

U_tot_mean = sum(U_tot)/len(U_tot)
print(U_tot_mean)

#%%

mode_energy = np.genfromtxt('energy_norm.csv', delimiter=',', skip_header=1)
timesteps = len(mode_energy[0,:])

fig, ax = plt.subplots()

for i in range(timesteps):    
    ax.plot(time, mode_energy[:, i])

ax.set_xlabel('time (arb.unit)')
ax.set_ylabel('signal (arb.unit)')  
ax.grid()
