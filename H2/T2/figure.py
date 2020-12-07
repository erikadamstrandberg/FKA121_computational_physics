#%% T! Benchmarking for variational Monte Carlo
# 
# Standard includes
import numpy as np
import matplotlib.pyplot as plt

#%% Uncorrelated! With only one electrion wavefunctions

array = np.genfromtxt('../T2/data/local_energy_bad_initial.csv', delimiter=',')

El = array
fig, ax = plt.subplots()
ax.plot(El, color='blue', label=r'$E_L[R_i]$')

ax.set_title(r'Local energy with bad initial step', fontsize='16')
ax.set_xlabel(r'$i$', fontsize='16')
ax.set_ylabel(r'$E_L$ [$a.u.$]', fontsize='16')

ax.grid()
ax.legend(fontsize='16', loc='upper right')

plt.savefig('../T2/figure/local_energy_bad_inital.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%%

El_avg = np.zeros(len(El)-1)
for i in range(len(El)-1):
    El_avg[i] = np.sum(El[0:i+1])/(i+1)


fig, ax = plt.subplots()
ax.plot(El_avg, color='orange', label=r'$E_L[R_i]$', linewidth=2)

ax.set_title(r'Calculated energy, $E_L$, with bad initial step', fontsize='16')
ax.set_xlabel(r'$i$', fontsize='16')
ax.set_ylabel(r'$E_L$ [$a.u.$]', fontsize='16')

ax.grid()
ax.legend(fontsize='16', loc='upper left')

plt.savefig('../T2/figure/energy_bad_inital.pdf', format='pdf', bbox_inches='tight')
plt.show()


#%%

array = np.genfromtxt('../T2/data/local_energy_bad_initial_burn_in.csv', delimiter=',')

El = array
fig, ax = plt.subplots()
ax.plot(El, color='blue', label=r'$E_L[R_i]$')

ax.set_title(r'Local energy with burn in, $N_{eq} = 1000$', fontsize='16')
ax.set_xlabel(r'$i$', fontsize='16')
ax.set_ylabel(r'$E_L$ [$a.u.$]', fontsize='16')

ax.grid()
ax.legend(fontsize='16', loc='upper right')

plt.savefig('../T2/figure/local_energy_bad_inital_burn_in.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%%

El_avg = np.zeros(len(El))
for i in range(len(El)):
    El_avg[i] = np.sum(El[0:i])/(i+1)


fig, ax = plt.subplots()
ax.plot(El_avg, color='orange', label=r'$E_L[R_i]$', linewidth=2)

ax.set_title(r'Calculated energy, $E_L$, with burn in $N_{eq} = 1000$', fontsize='16')
ax.set_xlabel(r'$i$', fontsize='16')
ax.set_ylabel(r'$E_L$ [$a.u.$]', fontsize='16')

ax.set_ylim([-3.2,-2.5])

ax.grid()
ax.legend(fontsize='16', loc='upper left')

plt.savefig('../T2/figure/energy_bad_inital_burn_in.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%%

array = np.genfromtxt('../T2/data/local_energy_good_initial.csv', delimiter=',')

El = array[0:4000]
fig, ax = plt.subplots()
ax.plot(El, color='blue', label=r'$E_L[R_i]$')

ax.set_title(r'Local energy with random initial and burn in, $N_{eq}=1000$', fontsize='16')
ax.set_xlabel(r'$i$', fontsize='16')
ax.set_ylabel(r'$E_L$ [$a.u.$]', fontsize='16')

ax.grid()
ax.legend(fontsize='16', loc='upper right')

plt.savefig('../T2/figure/local_energy_good_inital_burn_in.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%%

El_avg = np.zeros(len(El))
for i in range(len(El)):
    El_avg[i] = np.sum(El[0:i])/(i+1)


fig, ax = plt.subplots()
ax.plot(El_avg, color='orange', label=r'$E_L[R_i]$', linewidth=2)

ax.set_title(r'Calculated energy, $E_L$, with random initial and burn in $N_{eq} = 1000$', fontsize='16')
ax.set_xlabel(r'$i$', fontsize='16')
ax.set_ylabel(r'$E_L$ [$a.u.$]', fontsize='16')

ax.set_ylim([-3.2,-2.5])

ax.grid()
ax.legend(fontsize='16', loc='upper left')

plt.savefig('../T2/figure/energy_good_inital_burn_in.pdf', format='pdf', bbox_inches='tight')
plt.show()


#%%
array = np.genfromtxt('../T2/data/correlation.csv', delimiter=',')
N = len(array)
N = np.arange(0,N)

limit_for_ns = 0.135

x = np.array([np.min(N), np.max(N)])
y = np.array([limit_for_ns, limit_for_ns])

index_over = np.sum(array > limit_for_ns)
index_before = index_over-1

phi_over = array[index_over]
phi_before = array[index_before]
ns = index_over*(phi_before - limit_for_ns)/(phi_before-phi_over) + index_before*(limit_for_ns - phi_over)/(phi_before-phi_over)

x_interpol = np.array([ns, ns])
y_interpol = np.array([min(array), max(array)])

fig, ax = plt.subplots()
ax.plot(N, array, color='black', label=r'$\Phi_k$', linewidth=2)
ax.plot(x, y, '--', color='red', label=r'$\Phi_k = 0.135$', linewidth=2)
ax.plot(x_interpol, y_interpol, '--', label=r'$ns = 9.75$', color='blue', linewidth=2)

ax.set_title(r'Correlation function, $\Phi_k$', fontsize='16')
ax.set_xlabel(r'$k$', fontsize='16')
ax.set_ylabel(r'$\Phi_k$', fontsize='16')

ax.grid()
ax.legend(fontsize='16', loc='upper right')

plt.savefig('../T2/figure/ns_from_correlation.pdf', format='pdf', bbox_inches='tight')
plt.show()


#%%
array = np.genfromtxt('../T2/data/block_average.csv', delimiter=',')

N = len(array)
N = np.arange(1,N+1)

converged_b = array[100:-1]
mean = np.sum(converged_b)/(len(converged_b))
x = np.array([100, 500])
y = np.array([mean, mean])

fig, ax = plt.subplots()
ax.plot(N, array, color='black', label=r'$ns(B)$')
ax.plot(x,y, '--', color='red', linewidth=3, label=r'$ns = 9.76$')

ax.set_title(r'Block averaging', fontsize='16')
ax.set_xlabel(r'$B$', fontsize='16')
ax.set_ylabel(r'$ns$', fontsize='16')

ax.grid()
ax.legend(fontsize='16', loc='lower right')

plt.savefig('../T2/figure/ns_from_block.pdf', format='pdf', bbox_inches='tight')
plt.show()
