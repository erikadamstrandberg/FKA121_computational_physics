import numpy as np
import matplotlib.pyplot as plt

kb = 1.38064e-23
eV = 1.602e-19
N = 256

#%%

array = np.genfromtxt('data/energy_dt_1fs.csv', delimiter=',', skip_header=1)

scaling_const = 800

E_kin = array[:, 0]
E_pot = array[:, 1] + scaling_const
E_tot = array[:, 0] + array[:, 1] + scaling_const
time  = array[:, 2]


fig, ax = plt.subplots()
ax.plot(time, E_kin, color='red', label=r'$E_{kin}(t)$')
ax.plot(time, E_pot, color='blue', label=r'$E_{pot}(t)$')
ax.plot(time, E_tot, color='black', label=r'$E_{tot}(t)$')
    
ax.set_title(r'Conservation of energy, $dt=1$ $fs$ ', fontsize='16')
ax.set_xlabel(r'$t$ [$ps$]', fontsize='16')
ax.set_ylabel(r'$E$ [$arb. u.$]', fontsize='16')
ax.grid()
ax.legend(fontsize='16')

ax.set_xlim([0, 6])
ax.set_yticks([])

plt.savefig('figure/energy_dt_1fs.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%%
array = np.genfromtxt('data/energy_dt_10fs.csv', delimiter=',', skip_header=1)

E_kin = array[:, 0]
E_pot = array[:, 1] + scaling_const
E_tot = array[:, 0] + array[:, 1] + scaling_const
time  = array[:, 2] 


fig, ax = plt.subplots()
ax.plot(time, E_kin, color='red', label=r'$E_{kin}(t)$')
ax.plot(time, E_pot, color='blue', label=r'$E_{pot}(t)$')
ax.plot(time, E_tot, color='black', label=r'$E_{tot}(t)$')
    
ax.set_title(r'Conservation of energy, $dt=10$ $fs$', fontsize='16')
ax.set_xlabel(r'$t$ [$ps$]', fontsize='16')
ax.set_ylabel(r'$E$ [$arb. u.$]', fontsize='16')
ax.grid()
ax.legend(fontsize='16')

ax.set_xlim([0, 6])
ax.set_yticks([])

plt.savefig('figure/energy_dt_10fs.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%%
array = np.genfromtxt('data/energy_dt_20fs.csv', delimiter=',', skip_header=1)

E_kin = array[:, 0]
E_pot = array[:, 1] + scaling_const
E_tot = array[:, 0] + array[:, 1] + scaling_const
time  = array[:, 2]


fig, ax = plt.subplots()
ax.plot(time, E_kin, color='red', label=r'$E_{kin}(t)$', linewidth='4')
ax.plot(time, E_pot, color='blue', label=r'$E_{pot}(t)$')
ax.plot(time, E_tot, '--', color='black', label=r'$E_{tot}(t)$', linewidth='4')
    
ax.set_title(r'Conservation of energy, $dt=20$ $fs$', fontsize='16')
ax.set_xlabel(r'$t$ [$ps$]', fontsize='16')
ax.set_ylabel(r'$E$ [$arb. u.$]', fontsize='16')
ax.grid()
ax.legend(fontsize='16')

ax.set_xlim([0, 6])
ax.set_yticks([])

plt.savefig('figure/energy_dt_20fs.pdf', format='pdf', bbox_inches='tight')
plt.show()


#%%

array = array = np.genfromtxt('pos_after_verlet.csv', delimiter=',', skip_header=1)
x = array[:,0]
y = array[:,1]
z = array[:,2]

ax = plt.axes(projection='3d')
ax.scatter3D(x, y, z, color='black')

#%% 

array = array = np.genfromtxt('../T2/temp_and_pressure.csv', delimiter=',', skip_header=1)



