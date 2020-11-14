#!/usr/bin/env python
###############################################################################
# Create figure
# plt.savefig('name_of_figure.pdf', bbox_inches='tight')
###############################################################################

import numpy as np
import matplotlib.pyplot as plt


#%%

array = np.genfromtxt('energy.csv', delimiter=',', skip_header=1)
# U_kin = array[:, 0]
# U_pot = array[:, 1]
# time = array[:, 2]

fig, ax = plt.subplots()
ax.plot(array[:, 2], array[:, 0])
ax.plot(array[:, 2], array[:, 1])
ax.plot(array[:, 2], array[:, 0] + array[:,1])

ax.set_xlabel('time (arb.unit)')
ax.set_ylabel('signal (arb.unit)')
ax.grid()


