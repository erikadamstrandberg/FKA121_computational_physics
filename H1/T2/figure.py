import numpy as np
import matplotlib.pyplot as plt

kb = 1.306e-23
eV = 1.603e-19

#%%

N = 256

array = np.genfromtxt('energy.csv', delimiter=',', skip_header=1)
E_kin = array[:, 0]
E_pot = array[:, 1]
E_tot = array[:, 0] + array[:, 1]
time  = array[:, 2]


fig, ax = plt.subplots()
ax.plot(time, E_kin, color='red', label=r'$E_{kin}$')
ax.plot(time, E_pot, color='blue', label=r'$E_{pot}$')
ax.plot(time, E_tot, color='black', label=r'$E_{tot}$')

E_kin_mean = (sum(E_kin)/len(E_kin))*eV
T = (2/(3*N*kb))*E_kin_mean
print(T)


#%%
ax.set_title(r'$E_{pot}$ ', fontsize='22')
ax.set_xlabel(r'$V_{cell}$ [$Å^3$]', fontsize='16')
ax.set_ylabel(r'$E$ [$eV/unit\ cell$]', fontsize='16')
ax.grid()
ax.legend(fontsize='16')

#index_min = np.argmin(E_pot_fit)
#min_a = a_fit[index_min]

#plt.savefig('lattice_constant.pdf', format='pdf', bbox_inches='tight')
#plt.show()

#%%

array = np.genfromtxt('initial_random_displacement.csv', delimiter=',', skip_header=1)
x = array[:,0]
y = array[:,1]
z = array[:,2]

ax = plt.axes(projection='3d')
ax.scatter3D(x, y, z, color='black')

#%%

array = array = np.genfromtxt('after_verlet.csv', delimiter=',', skip_header=1)
x = array[:,0]
y = array[:,1]
z = array[:,2]

ax = plt.axes(projection='3d')
ax.scatter3D(x, y, z, color='black')
