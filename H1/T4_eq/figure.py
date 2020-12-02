import numpy as np
import matplotlib.pyplot as plt

kb = 1.38064e-23
eV = 1.602e-19
N = 256

#%%

array = np.genfromtxt('../T4_eq/data/TPV_700C_1fs_1.csv', delimiter=',', skip_header=1)

temp     = array[:, 0] - 272.15
pressure = array[:, 1]
volume   = array[:, 2]
time     = array[:, 3]

#%% Plotting temp equil

T_equil = 700
T_melt = 1600
T_start = 2
T_start_2 = 16
T_start_3 = 24

y_min = 350
y_max = 1800

fig, ax = plt.subplots()
ax.plot(time, temp, color='blue', label=r"$T(t)$")
ax.plot([min(time), max(time)], [T_equil, T_equil], color='black', label=r"$T_{liquid}=700^{\circ}$C", linewidth='4')
ax.plot([min(time), max(time)], [T_melt, T_melt], color='black', label=r"$T_{melt}=1600^{\circ}$C", linewidth='2')
ax.plot([T_start, T_start], [y_min, y_max], '--',   color='red', label=r"$t_{T,1}= 2\ ps$", linewidth='2')
ax.plot([T_start_2, T_start_2], [y_min, y_max], '--',       color= 'red' , label=r"$t_{T,2}= 16\ ps$", linewidth='2')
ax.plot([T_start_3, T_start_3], [y_min, y_max], '--',       color= 'red' , label=r"$t_{T,3}= 24\ ps$", linewidth='2')

ax.set_title(r'Equilibration of temperature, $T_{liquid}=700\ ^{\circ}C$', fontsize='16')
ax.set_xlabel(r'$t$ [$ps$]', fontsize='16')
ax.set_ylabel(r'$T$ [$^{\circ}C$]', fontsize='16')
ax.grid()
ax.legend(fontsize='12')

ax.set_xlim([0, 40])
ax.set_ylim([y_min, y_max])

plt.savefig('figure/temp_equil_T700.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%% Plotting pressure equil

P_equil = 1e-4
T_start = 8
T_start_2 = 18
T_start_3 = 26

y_min = -1
y_max = 10

fig, ax = plt.subplots()
ax.plot(time, pressure, color='orange', label=r"$P(t)$")
ax.plot([min(time), max(time)], [P_equil, P_equil], color='black', label=r"$P_{liquid}$")
ax.plot([T_start, T_start], [y_min, y_max], '--',   color='red', label=r"$t_{P,1}= 8$ ps", linewidth='2')
ax.plot([T_start_2,T_start_2], [y_min, y_max], '--',       color= 'red' , label=r"$t_{P,2}= 18$ ps", linewidth='2')
ax.plot([T_start_3,T_start_3], [y_min, y_max], '--',       color= 'red' , label=r"$t_{P,3}= 26$ ps", linewidth='2')

ax.set_title(r'Equilibration of pressure, $P_{liquid}=10^{-4}$ GPa ', fontsize='16')
ax.set_xlabel(r'$t$ [$ps$]', fontsize='16')
ax.set_ylabel(r'$P$ [$GPa$]', fontsize='16')
ax.grid()
ax.legend(fontsize='16')
ax.set_ylim([y_min, y_max])
ax.set_xlim([0, 40])

plt.savefig('figure/pressure_equil_T700.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%% Volume 

fig, ax = plt.subplots()
ax.plot(time, volume, color='black', label=r'$V(t)$')

ax.set_title(r'Volume during equilibrium for solid phase ', fontsize='16')
ax.set_xlabel(r'$t$ [$ps$]', fontsize='16')
ax.set_ylabel(r'$V$ [$Å^3$]', fontsize='16')

ax.grid()
ax.legend(fontsize='16')
ax.set_xlim([0, 40])

plt.savefig('figure/volume_equil_T700.pdf', format='pdf', bbox_inches='tight')
plt.show()
