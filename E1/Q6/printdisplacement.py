#!/usr/bin/env python
###############################################################################
# E1code2
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv
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

plt.show()

#%%

array = np.genfromtxt('powerspectrum.csv', delimiter=',', skip_header=1)
omega0 = np.genfromtxt('omega0.csv', delimiter=',', skip_header=1)

power_spectrum = array[:, 0]+array[:, 1]+array[:, 2]
scale = omega0/(2*np.pi)

x_min = -100/scale
x_max = 100/scale

y_min = -0.01
y_max = max(power_spectrum)

fig, ax = plt.subplots()
ax.plot(array[:, 3]/scale, power_spectrum)

ax.set_xlabel(r'$\omega_0$')
ax.set_ylabel(r'power spectrum [arb.u]')
ax.grid()

x_lim = np.array([x_min,x_max])
y_lim = np.array([y_min,y_max])
ax.set_xlim(x_lim)
ax.set_ylim(y_lim)