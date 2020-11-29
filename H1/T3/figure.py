import numpy as np
import matplotlib.pyplot as plt

kb = 1.38064e-23
eV = 1.602e-19
N = 256

#%%

array = np.genfromtxt('energy.csv', delimiter=',', skip_header=1)

E_kin = array[:, 0]
E_pot = array[:, 1]
E_tot = E_kin + E_pot
time = array[:, 2]

fig, ax = plt.subplots()
ax.plot(time, E_kin)
ax.plot(time, E_pot)
ax.plot(time, E_tot)

#%%

array = np.genfromtxt('TPV.csv', delimiter=',', skip_header=1)

start = 1
N = 1000

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
    
    
x = np.array([0,100])
y = np.array([0.0001,0.0001])

fig, ax = plt.subplots()
ax.plot(time, P_average)
ax.plot(x, y)

ax.set_xlim([0, 100])
ax.set_ylim([-0.0001, 0.02])