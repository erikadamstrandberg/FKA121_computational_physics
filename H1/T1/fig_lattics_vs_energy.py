import numpy as np
import matplotlib.pyplot as plt


#%%

array = np.genfromtxt('lattice_pot_energy.csv', delimiter=',', skip_header=1)
a = array[:, 0]
E_pot = array[:, 1]

fit_coeff = np.polyfit(a, E_pot, 2)
a_fit = np.linspace(a[0], a[-1], 1000);
E_pot_fit = fit_coeff[0]*a_fit**2 + fit_coeff[1]*a_fit + fit_coeff[2]


fig, ax = plt.subplots()
ax.plot(a**3, E_pot, 'x', color='black', label=r'$E_{sim}$')
ax.plot(a_fit**3, E_pot_fit, color='red', linewidth='0.8', label=r'$E_{fit}$')

ax.set_title(r'$E_{pot}$ ', fontsize='16')
ax.set_xlabel(r'$V_{cell}$ [$Å^3$]', fontsize='16')
ax.set_ylabel(r'$E$ [$eV/unit\ cell$]', fontsize='16')
ax.grid()
ax.legend(fontsize='16')

index_min = np.argmin(E_pot_fit)
min_a = a_fit[index_min]

plt.savefig('lattice_constant.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%%

array = array = np.genfromtxt('initial_AL.csv', delimiter=',', skip_header=1)
x = array[:,0]
y = array[:,1]
z = array[:,2]

ax = plt.axes(projection='3d')
ax.scatter3D(x, y, z, color='black')

