import numpy as np
import matplotlib.pyplot as plt

kb = 1.38064e-23
eV = 1.602e-19
N = 256

#%%

array = np.genfromtxt('energy.csv', delimiter=',', skip_header=1)

E_kin = array[:, 0]
E_pot = array[:, 1]
E_tot = array[:, 0] + array[:, 1]
time  = array[:, 2]

fig, ax = plt.subplots()
ax.plot(time, E_kin)
ax.plot(time, E_pot)
ax.plot(time, E_tot)


#%%

array = np.genfromtxt('temp_and_pressure.csv', delimiter=',', skip_header=1)

temp     = array[:, 0]
pressure = array[:, 1]
time     = array[:, 2]

fig, ax = plt.subplots()
ax.plot(time, temp*1e-2)
ax.plot(time, pressure*160.2)


#%%
T_equil = 500e-2
ax.plot([min(time), max(time)], [T_equil, T_equil])

P_equil = 1e-4
ax.plot([min(time), max(time)], [P_equil, P_equil])

tau_p = 150*1e-3 + 0.5
tau_a = 200*1e-3 + 0.5
ax.plot([tau_p, tau_p], [0,10])
ax.plot([tau_a, tau_a], [0,10])

ax.set_xlim([0.5,5])
ax.set_ylim([-0.02,1e4])
    
#%%

array = array = np.genfromtxt('initial_random_displacement.csv', delimiter=',', skip_header=1)
x = array[:,0]
y = array[:,1]
z = array[:,2]

ax = plt.axes(projection='3d')
ax.scatter3D(x, y, z, color='black')
