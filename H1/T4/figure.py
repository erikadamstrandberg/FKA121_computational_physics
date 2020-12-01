import numpy as np
import matplotlib.pyplot as plt

kb = 1.38064e-23
eV = 1.602e-19
N = 256

#%%

array = np.genfromtxt('../T4/energy.csv', delimiter=',', skip_header=1)

E_kin = array[:, 0]
E_pot = array[:, 1]
E_tot = E_kin + E_pot
time = array[:, 2]

fig, ax = plt.subplots()
ax.plot(time, E_kin)
ax.plot(time, E_pot)
ax.plot(time, E_tot)

#%%

array = np.genfromtxt('../T4/TPV.csv', delimiter=',', skip_header=1)

start = 1
N = len(array)

temp     = array[start:N, 0]
pressure = array[start:N, 1]
time     = array[start:N, 2]

T_average = np.zeros(len(temp))
for i in range(len(T_average)):
    T_average[i] = np.sum(temp[0:i])/i
 
fig, ax = plt.subplots()
ax.plot(time, T_average)

P_average = np.zeros(len(temp))
for i in range(len(P_average)):
    P_average[i] = np.sum(pressure[0:i])/i
    
    
x = np.array([0,20])
y = np.array([0.0001,0])
    
fig, ax = plt.subplots()
ax.plot(time, P_average)

#ax.set_xlim([0, 100])
#ax.set_ylim([-0.01, 0.02])


#%% Timetrails

pos_data = np.genfromtxt('../T4/q_trail.csv', delimiter=',')
M, n_particles = np.shape(pos_data)
n_particles = int(n_particles/3)

atom_number = 2

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(pos_data[:,atom_number], pos_data[:,atom_number+1], pos_data[:,atom_number+2])



