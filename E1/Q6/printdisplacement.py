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


for i in range(len(omega)):
    x = np.array([omega[i], omega[i]])/(2*np.pi)
    y = np.array([y_min, y_max])
    #ax.plot(x,y)
    
    
fund_cen_1 = 2349.0*1e2
fund_wavelength_1 = 1/fund_cen_1
fund_freq_1 = (3e8/fund_wavelength_1)*1e-12

fund_cen_2 = 1333*1e2
fund_wavelength_2 = 1/fund_cen_2
fund_freq_2 = (3e8/fund_wavelength_2)*1e-12

x
freq_measured = np.array([fund_freq_1, fund_freq_2])

for i in range(len(freq_measured)):
    x = np.array([freq_measured[i], freq_measured[i]])
    y = np.array([y_min, y_max])
    ax.plot(x,y)
    
plt.savefig('powerspectrum_with_measured.pdf', bbox_inches='tight')


#%%