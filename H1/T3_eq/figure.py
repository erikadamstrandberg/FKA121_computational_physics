import numpy as np
import matplotlib.pyplot as plt

kb = 1.38064e-23
eV = 1.602e-19
N = 256

#%%

array = np.genfromtxt('../T3_eq/data/TPV_500C_1fs_1.csv', delimiter=',', skip_header=1)

temp     = array[:, 0] - 272.15
pressure = array[:, 1]
volume   = array[:, 2]
time     = array[:, 3]

#%% Plotting temp equil

T_equil = 500
T_start = 2
T_start_2 = 8

y_min = 600 - 272.15
y_max = 1100 - 272.15

fig, ax = plt.subplots()
ax.plot(time, temp, color='blue', label=r"$T(t)$")
ax.plot([min(time), max(time)], [T_equil, T_equil], color='black', label=r"$T_{solid}$", linewidth='4')
ax.plot([T_start, T_start], [y_min, y_max], '--',   color='red', label=r"$t_{T,1}= 2$ $ps$",linewidth='2')
ax.plot([T_start_2, T_start_2], [y_min, y_max], '--',   color='red', label=r"$t_{T,2}= 8$ $ps$",linewidth='2')

ax.set_title(r'Equilibration of temperature, $T_{solid}=500\ ^{\circ}C$', fontsize='16')
ax.set_xlabel(r'$t$ [$ps$]', fontsize='16')
ax.set_ylabel(r'$T$ [$^{\circ}C$]', fontsize='16')
ax.grid()
ax.legend(fontsize='16', loc='upper right')

ax.set_xlim([0, 20])
ax.set_ylim([y_min, y_max])

plt.savefig('figure/temp_equil_T500.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%% Plotting pressure equil

P_equil = 1e-4
T_start = 6
T_start_2 = 10

y_min = -0.3
y_max = 5

fig, ax = plt.subplots()
ax.plot(time, pressure, color='orange', label=r"$P(t)$")
ax.plot([min(time), max(time)], [P_equil, P_equil], color='black', label=r"$P_{solid}$")
ax.plot([T_start, T_start], [y_min, y_max], '--',   color='red', label=r"$t_{P,1}= 6$ $ps$")
ax.plot([T_start_2, T_start_2], [y_min, y_max], '--',       color= 'red' , label=r"$t_{p,2} = 10\ ps$")

ax.set_title(r'Equilibration of pressure, $P_{solid}=10^{-4}$ $GPa$', fontsize='16')
ax.set_xlabel(r'$t$ [$ps$]', fontsize='16')
ax.set_ylabel(r'$P$ [$GPa$]', fontsize='16')
ax.grid()
ax.legend(fontsize='16')

ax.set_xlim([0, 20])
ax.set_ylim([y_min, y_max])

plt.savefig('figure/pressure_equil_T500.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%% Volume 

fig, ax = plt.subplots()
ax.plot(time, volume, color='black', label=r'$V(t)$')

ax.set_title(r'Volume during equilibrium for solid phase ', fontsize='16')
ax.set_xlabel(r'$t$ [$ps$]', fontsize='16')
ax.set_ylabel(r'$V$ [$Å^3$]', fontsize='16')

ax.grid()
ax.legend(fontsize='16')

ax.set_xlim([0, 20])

plt.savefig('figure/volume_equil_T500.pdf', format='pdf', bbox_inches='tight')
plt.show()

