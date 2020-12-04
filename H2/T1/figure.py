#%% T! Benchmarking for variational Monte Carlo
# 
# Standard includes

import numpy as np
import matplotlib.pyplot as plt


#%% Centeral mean field approx. from Hartree

Z_unscreened = 2
Z_optimized  = 27/16

r = np.linspace(0, 5, 1000)
def rho(r, Z):
    return Z**3*4*r**2*np.exp(-2*Z*r)

# NORMALIZE!?!
rho_unscreened = rho(r, Z_unscreened)
rho_optimized = rho(r, Z_optimized)
fig, ax = plt.subplots()
ax.plot(r, rho_unscreened, '--', color='blue', label=r'$\rho_{unscreend}$')
ax.plot(r, rho_optimized , '--', color='red', label=r'$\rho_{optimized}$')

ax.set_title(r'Probability of finding a electron', fontsize='16')
ax.set_xlabel(r'$r$ [$Å$]', fontsize='16')
ax.set_ylabel(r'$\rho$', fontsize='16')

ax.grid()
ax.legend(fontsize='16', loc='upper right')

plt.savefig('../T1/figure/probability_benchmark.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%%


