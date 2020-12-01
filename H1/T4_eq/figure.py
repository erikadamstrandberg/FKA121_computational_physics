import numpy as np
import matplotlib.pyplot as plt

kb = 1.38064e-23
eV = 1.602e-19
N = 256

#%%

array = np.genfromtxt('../T4_eq/data/TPV_700C_1fs_1.csv', delimiter=',', skip_header=1)

temp     = array[:, 0]
pressure = array[:, 1]
volume   = array[:, 2]
time     = array[:, 3]

#%% Plotting temp equil

T_equil = 700 + 272.15
T_start = 2
T_tau   = T_start + 200*1e-3

y_min = 650
y_max = 1500

fig, ax = plt.subplots()
ax.plot(time, temp, color='blue', label=r"$T(t)$")
ax.plot([min(time), max(time)], [T_equil, T_equil], color='black', label=r"$T_{solid}=700^{\circ}$C")
ax.plot([T_start, T_start], [y_min, y_max], '--',   color='black', label=r"$t_{start}= 2$ ps",)
ax.plot([T_tau, T_tau], [y_min, y_max], '--',       color= 'red' , label=r"$\tau_{t}= 200$ dt")

ax.set_title(r'Equilibration of temperature, $T_{solid}=500$ K', fontsize='16')
ax.set_xlabel(r'$t$ [ps]', fontsize='16')
ax.set_ylabel(r'$T$ [K]', fontsize='16')
ax.grid()
ax.legend(fontsize='16')

ax.set_ylim([y_min, y_max])

plt.savefig('figure/temp_equil_T700.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%% Plotting pressure equil

P_equil = 1e-4
T_start = 4
T_tau = T_start + 400*3e-3

y_min = -1
y_max = 5

fig, ax = plt.subplots()
ax.plot(time, pressure, color='orange', label=r"$P(t)$")
ax.plot([min(time), max(time)], [P_equil, P_equil], color='black', label=r"$P_{solid}$")
ax.plot([T_start, T_start], [y_min, y_max], '--',   color='black', label=r"$t_{start}= 4$ ps")
ax.plot([T_tau, T_tau], [y_min, y_max], '--',       color= 'red' , label=r"$\tau_{p}= 400$ dt")

ax.set_title(r'Equilibration of pressure, $P_{solid}=10^{-4}$ GPa ', fontsize='16')
ax.set_xlabel(r'$t$ [ps]', fontsize='16')
ax.set_ylabel(r'$P$ [GPa]', fontsize='16')
ax.grid()
ax.legend(fontsize='16')
ax.set_ylim([y_min, y_max])

plt.savefig('figure/pressure_equil_T700.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%% Volume 

fig, ax = plt.subplots()
ax.plot(time, volume, color='black', label=r'$V(t)$')

ax.set_title(r'Volume during equilibrium for solid phase ', fontsize='16')
ax.set_xlabel(r'$t$ [ps]', fontsize='16')
ax.set_ylabel(r'$V$ [$Å^3$]', fontsize='16')

ax.grid()
ax.legend(fontsize='16')

plt.savefig('figure/volume_equil_T700.pdf', format='pdf', bbox_inches='tight')
plt.show()
