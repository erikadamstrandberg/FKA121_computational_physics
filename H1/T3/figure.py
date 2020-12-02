import numpy as np
import matplotlib.pyplot as plt

kb = 1.38064e-23
eV = 1.602e-19
N = 256


#%%

array = np.genfromtxt('../T3/data/merged_TPV.csv', delimiter=',')

temp     = array[:, 0] - 272.15
pressure = array[:, 1]
time     = array[:, 2]

dt = time[1] - time[0]
T = len(time)*dt
time = np.arange(0,T,dt)


T_average = np.zeros(len(temp))
for i in range(len(T_average)):
    T_average[i] = np.sum(temp[0:i])/i
 
    
x = np.array([0, 5])
y = np.array([500, 500])

fig, ax = plt.subplots()
ax.plot(time-1, T_average, color='blue', label=r"$T_{avg}(t)$")
ax.plot(x,y, '--', color='black', label=r"$T_{solid}$")

ax.set_title(r'Average temperature after equilibration, $T_{solid}=500\ ^{\circ}C$ ', fontsize='16')
ax.set_xlabel(r'$t$ [$ps$]', fontsize='16')
ax.set_ylabel(r'$T_{avg}$ [$^{\circ}C$]', fontsize='16')
ax.grid()
ax.legend(fontsize='16', loc='lower right')

ax.set_xlim([0, 4])
ax.set_ylim([460, 510])

plt.savefig('figure/temp_equil_T500_prod.pdf', format='pdf', bbox_inches='tight')
plt.show()


#%%

P_average = np.zeros(len(temp))
for i in range(len(P_average)):
    P_average[i] = np.sum(pressure[0:i])/i
    
    
x = np.array([0,100])
y = np.array([0.0001,0.0001])

fig, ax = plt.subplots()
ax.plot(time-1, P_average, color='orange', label=r'$P_{avg}(t)$')

ax.set_title(r'Average pressure after equilibration, $P_{solid}=10^{-4}\ GPa$ ', fontsize='16')
ax.set_xlabel(r'$t$ [$ps$]', fontsize='16')
ax.set_ylabel(r'$P_{avg}$ [$GPa$]', fontsize='16')
ax.grid()
ax.legend(fontsize='16', loc='lower right')

ax.set_xlim([0, 4])
ax.set_ylim([0, 0.1])

plt.savefig('figure/pressure_equil_T500_prod.pdf', format='pdf', bbox_inches='tight')
plt.show()


#%% Timetrails

pos_data = np.genfromtxt('../T3/data/merged_q_trail.csv', delimiter=',')
M, n_particles = np.shape(pos_data)
n_particles = int(n_particles/3)

atom_number = 1

x = pos_data[:,atom_number]
y = pos_data[:,atom_number+1]
z = pos_data[:,atom_number+2]

x = x - x[0]
y = y - y[0]
z = z - z[0]

fig, ax = plt.subplots()
ax.plot(time, x, color='black', label=r'$x(t)$')
ax.set_title(r'Displacement of atom in solid phase', fontsize='16')
ax.set_xlabel(r'$t$ [$ps$]', fontsize='16')
ax.set_ylabel(r'$x$ [$Å$]', fontsize='16')
ax.grid()
ax.legend(fontsize='16', loc='lower right')

ax.set_xlim([0, 4])

plt.savefig('figure/x_displacement_T500.pdf', format='pdf', bbox_inches='tight')
plt.show()

fig, ax = plt.subplots()
ax.plot(time, y, color='blue', label=r'$y(t)$')
ax.set_title(r'Displacement of atom in solid phase', fontsize='16')
ax.set_xlabel(r'$t$ [$ps$]', fontsize='16')
ax.set_ylabel(r'$y$ [$Å$]', fontsize='16')
ax.grid()
ax.legend(fontsize='16', loc='lower right')

ax.set_xlim([0, 4])

plt.savefig('figure/y_displacement_T500.pdf', format='pdf', bbox_inches='tight')
plt.show()

fig, ax = plt.subplots()
ax.plot(time, z, color='red', label=r'$z(t)$')
ax.set_title(r'Displacement of atom in solid phase', fontsize='16')
ax.set_xlabel(r'$t$ [$ps$]', fontsize='16')
ax.set_ylabel(r'$z$ [$Å$]', fontsize='16')
ax.grid()
ax.legend(fontsize='16', loc='lower right')

ax.set_xlim([0, 4])

plt.savefig('figure/z_displacement_T500.pdf', format='pdf', bbox_inches='tight')
plt.show()