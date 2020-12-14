#%% T! Benchmarking for variational Monte Carlo
# 
# Standard includes
import numpy as np
import matplotlib.pyplot as plt

#%%
def read_data(alpha, what_data):
    data = np.genfromtxt('../T3_long/data/finalvalue_{:.2f}.csv'.format(alpha), delimiter=',')
#   E_l = array[0]
#   sigma2 = array[1]
#   alpha = array[2]
#   N = array[3]
#   ns = array[4]
    return data[what_data]

def get_E_l(alpha):
    return read_data(alpha, 0)

def get_sigma2(alpha):
    return read_data(alpha, 1)

def get_alpha(alpha):
    return read_data(alpha, 2)

def get_N(alpha):
    return read_data(alpha, 3)

def get_ns(alpha):
    return read_data(alpha, 4)


a = np.array([0.05, 0.08, 0.10, 0.12, 0.13, 0.14, 0.15, 0.16, 0.18, 0.2, 0.22, 0.25])
Na     = len(a)
E_l    = np.zeros(Na)
sigma2 = np.zeros(Na)
alpha  = np.zeros(Na)
N      = np.zeros(Na)
ns     = np.zeros(Na)
error  = np.zeros(Na)

for i in range(Na):
    E_l[i] = get_E_l(a[i])
    sigma2[i] = get_sigma2(a[i])
    alpha[i] = get_alpha(a[i])
    N[i] = get_N(a[i])
    ns[i] = get_ns(a[i])
    error = np.sqrt(sigma2[i]*ns[i]/N[i])

fig, ax = plt.subplots()
ax.errorbar(alpha, E_l, yerr=error, fmt='--o', color='black', ecolor='black', label=r'$E(\alpha)$')



ax.set_title(r'Energy sweeping over alpha', fontsize='16')
ax.set_xlabel(r'$\alpha$', fontsize='16')
ax.set_ylabel(r'$E$ [$a.u.$]', fontsize='16')

ax.grid()
ax.legend(fontsize='16', loc='upper right')

plt.savefig('../T3_long/figure/alpha_sweep.pdf', format='pdf', bbox_inches='tight')
plt.show()

