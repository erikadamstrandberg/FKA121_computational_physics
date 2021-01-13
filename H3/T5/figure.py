import numpy as np
import matplotlib.pyplot as plt

c = 0.05

dir_name = 'data_' + str(c)

E_i = np.genfromtxt(dir_name + '/E_i.txt')

T_1 = np.genfromtxt(dir_name + '/T1.txt')
T_2 = np.genfromtxt(dir_name + '/T2.txt')
R_1 = np.genfromtxt(dir_name + '/R1.txt')
R_2 = np.genfromtxt(dir_name + '/R2.txt')

fig, ax = plt.subplots()
ax.plot(E_i, T_1, color='b', label=r'$T_1$')
ax.plot(E_i, R_1, '--', color='b', label=r'$R_1$')
ax.plot(E_i, T_2, color='r', label=r'$T_2$')
#ax.plot(E_i, R_2, label=r'$R_2$')
ax.set_xlabel(r'$E_{initial}$ [eV]', fontsize=16)
ax.set_ylabel(r'$T, R$', fontsize=16)
ax.set_title('Transmission and reflection coefficients, c = ' + str(c) + ' [eV]', fontsize=16)
ax.legend(fontsize=16)
ax.grid()

#plt.savefig('T_R_'+str(c)+'.pdf', format='pdf', bbox_inches='tight')

plt.show()

# plot probability densities

index_R1 = 0
index_T2 = 1
index_T1 = 2
index = index_T1

E_i = np.genfromtxt(dir_name + '/E_i_new.txt')
x = np.genfromtxt(dir_name + '/x.txt')

n_x_1 = np.genfromtxt(dir_name + '/n_x_1_' + str(index) + '.txt')
n_x_2 = np.genfromtxt(dir_name + '/n_x_2_' + str(index) + '.txt')

fig, ax = plt.subplots()
ax.plot(x, n_x_1, color='b', label=r'$|\psi_1(x)|^2$')
ax.plot(x, n_x_2, color='r', label=r'$|\psi_2(x)|^2$')
ax.set_xlabel(r'$x$ [Å]', fontsize=16)
ax.set_ylabel(r'$|\psi(x)|^2$ $[Å^{-1}]$', fontsize=16)
ax.set_title(r'Probability density, $E_{initial}$ = ' + str(E_i[index]) + ' [eV]', fontsize=16)
ax.legend(fontsize=16)
ax.grid()

#plt.savefig('n_x_'+str(E_i[index])+'.pdf', format='pdf', bbox_inches='tight')

plt.show()
