import numpy as np
import matplotlib.pyplot as plt

kb = 1.38064e-23
eV = 1.602e-19
N = 256

#%%

array = np.genfromtxt('data/TPV.csv', delimiter=',', skip_header=1)

temp     = array[:, 0]
pressure = array[:, 1]
volume   = array[:, 2]
time     = array[:, 3]

#%% Plotting temp equil

T_equil = 500
T_start = 2
T_tau   = T_start + 200*1e-3

y_min = 300
y_max = 1000

fig, ax = plt.subplots()
ax.plot(time, temp, color='blue', label=r"$T(t)$")
ax.plot([min(time), max(time)], [T_equil, T_equil], color='black', label=r"$T_{solid}=500$ K")
ax.plot([T_start, T_start], [y_min, y_max], '--',   color='black', label=r"$t_{start}= 2$ ps",)
ax.plot([T_tau, T_tau], [y_min, y_max], '--',       color= 'red' , label=r"$\tau_{t}= 200$ dt")

ax.set_title(r'Equilibration of temperature, $T_{solid}=500$ K', fontsize='16')
ax.set_xlabel(r'$t$ [ps]', fontsize='16')
ax.set_ylabel(r'$T$ [K]', fontsize='16')
ax.grid()
ax.legend(fontsize='16')

ax.set_ylim([y_min, y_max])

plt.savefig('figure/temp_equil_T500.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%% Plotting pressure equil

P_equil = 1e-4
T_start = 4
T_tau = T_start + 400*1e-3

y_min = -1
y_max = 5

pressure_GPa = pressure*160.2

fig, ax = plt.subplots()
ax.plot(time, pressure_GPa, color='orange', label=r"$P(t)$")
ax.plot([min(time), max(time)], [P_equil, P_equil], color='black', label=r"$P_{solid}$")
ax.plot([T_start, T_start], [y_min, y_max], '--',   color='black', label=r"$t_{start}= 4$ ps")
ax.plot([T_tau, T_tau], [y_min, y_max], '--',       color= 'red' , label=r"$\tau_{p}= 400$ dt")

ax.set_title(r'Equilibration of pressure, $P_{solid}=10^{-4}$ GPa ', fontsize='16')
ax.set_xlabel(r'$t$ [ps]', fontsize='16')
ax.set_ylabel(r'$P$ [GPa]', fontsize='16')
ax.grid()
ax.legend(fontsize='16')
ax.set_ylim([y_min, y_max])

plt.savefig('figure/pressure_equil_T500.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%% Volume 

fig, ax = plt.subplots()
ax.plot(time, volume, color='black', label=r'$V(t)$')

ax.set_title(r'Volume during equilibrium for solid phase ', fontsize='16')
ax.set_xlabel(r'$t$ [ps]', fontsize='16')
ax.set_ylabel(r'$V$ [$Å^3$]', fontsize='16')

ax.grid()
ax.legend(fontsize='16')

plt.savefig('figure/volume_equil_T500.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%% Timetrails

array = np.genfromtxt('timetrails.csv', delimiter=',', skip_header=1)
length_saved = len(array)
NDIM = 3
number_of_atoms = 5

q = np.zeros((length_saved, NDIM))

plot_from = 1000
plot_to = 4000

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for i in range(1):
    q[:,0] = array[:,0 + 3*i]
    q[:,1] = array[:,1 + 3*i]
    q[:,2] = array[:,2 + 3*i]

    ax.plot(q[plot_from:plot_to,0], q[plot_from:plot_to,1], q[plot_from:plot_to,2])
    
    

#%%
    
asasd = 0













