import numpy as np
import matplotlib.pyplot as plt


#%%

array = np.genfromtxt('lattice_pot_energy.csv', delimiter=',', skip_header=1)
a = array[:, 0]
E_pot = array[:, 1]


fig, ax = plt.subplots()
ax.plot(a**3, E_pot, 'x', color='black', label=r'$E_{sim}$')


ax.set_title(r'$E_{pot}$ ', fontsize='22')
ax.set_xlabel(r'$V_{cell}$ [$Å^3$]', fontsize='16')
ax.set_ylabel(r'$E$ [$eV/unit\ cell$]', fontsize='16')
ax.grid()
ax.legend(fontsize='16')

#index_min = np.argmin(E_pot_fit)
#min_a = a_fit[index_min]

#plt.savefig('lattice_constant.pdf', format='pdf', bbox_inches='tight')
#plt.show()

