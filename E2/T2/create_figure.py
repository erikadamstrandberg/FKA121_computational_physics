#!/usr/bin/env python
###############################################################################
# Create figure
# plt.savefig('name_of_figure.pdf', bbox_inches='tight')
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

#%%
array = np.genfromtxt('timetrail.csv', delimiter=',', skip_header=1)
# q1 = array[:, 0]
# q2 = array[:, 1]
# q3 = array[:, 2]
# U_kin = array[:, 3]
# U_pot = array[:, 4]
# time = array[:, 5]

fig, ax = plt.subplots()
ax.plot(array[:, 5], array[:, 0])
ax.plot(array[:, 5], array[:, 1])
ax.plot(array[:, 5], array[:, 2])

ax.set_xlabel('time (arb.unit)')
ax.set_ylabel('signal (arb.unit)')
ax.grid()

plt.savefig('displacement.pdf', bbox_inches='tight')

#%%

fig, ax = plt.subplots()
ax.plot(array[:, 5], array[:, 3])
ax.plot(array[:, 5], array[:, 4])
ax.plot(array[:, 5], array[:, 3] + array[:, 4])

ax.set_xlabel('time (arb.unit)')
ax.set_ylabel('power (arb.unit)')
ax.grid()

plt.savefig('energy.pdf', bbox_inches='tight')

#%%

array = np.genfromtxt('powerspectrum.csv', delimiter=',', skip_header=1)
omega0 = np.genfromtxt('omega0.csv', delimiter=',', skip_header=1)
omega = np.array([0,omega0[0],omega0[1]])

power_spectrum = array[:, 0]+array[:, 1]+array[:, 2]

fig, ax = plt.subplots()
ax.plot(array[:, 3], power_spectrum)

ax.set_xlabel(r'freq [THz]')
ax.set_ylabel(r'power spectrum [arb.u]')
ax.grid()

x_min = 74
x_max = 75

y_min = -0.01
y_max = max(power_spectrum)

x_lim = np.array([x_min,x_max])
y_lim = np.array([y_min,y_max])
ax.set_xlim(x_lim)
ax.set_ylim(y_lim)
