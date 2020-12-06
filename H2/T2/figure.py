#%% T! Benchmarking for variational Monte Carlo
# 
# Standard includes
import numpy as np
import matplotlib.pyplot as plt

#%% Uncorrelated! With only one electrion wavefunctions
array = np.genfromtxt('../T1/data/markov_chain_uncorrelated.csv', delimiter=',')

## Markox chain with coordinates for first electron
x1 = array[:,0]
y1 = array[:,1]
z1 = array[:,2]

## Markox chain with coordinates for second electron    
x2 = array[:,3]
y2 = array[:,4]
z2 = array[:,5]

r1 = np.zeros(len(array[:,0]))
r2 = np.zeros(len(array[:,0]))
theta = np.zeros(len(array[:,0]))
uniform = np.zeros(len(array[:,0]))

for i in range(len(r1)):
    r1[i] = np.sqrt(np.dot(x1[i],x1[i]) + np.dot(y1[i],y1[i]) + np.dot(z1[i],z1[i]))
    r2[i] = np.sqrt(np.dot(x2[i],x2[i]) + np.dot(y2[i],y2[i]) + np.dot(z2[i],z2[i]))
    
    arg = (x1[i]*x2[i] + y1[i]*y2[i] + z1[i]*z2[i])
    theta[i] = np.arccos(arg/(r1[i]*r2[i]))
    uniform[i] = np.cos(theta[i])


#%% Bins for histogram
step_bin = 0.1
start_bin = -1
stopp_bin = 1 + step_bin
b = np.arange(start_bin, stopp_bin, step_bin)

fig, ax = plt.subplots()
ax.hist(uniform, bins=b, density=True)


#%% Bins for histogram
step_bin = 0.1
start_bin = 0
stopp_bin = np.pi + step_bin
b = np.arange(start_bin, stopp_bin, step_bin)

fig, ax = plt.subplots()
ax.hist(theta, bins=b, density=True)


#%% Correlated data
array = np.genfromtxt('../T1/data/markov_chain_correlated.csv', delimiter=',')

## Markox chain with coordinates for first electron
x1 = array[:,0]
y1 = array[:,1]
z1 = array[:,2]

## Markox chain with coordinates for second electron    
x2 = array[:,3]
y2 = array[:,4]
z2 = array[:,5]

El = array[:,6]

N = len(array[:,0])


r1 = np.zeros(len(array[:,0]))
r2 = np.zeros(len(array[:,0]))
theta = np.zeros(len(array[:,0]))
uniform = np.zeros(len(array[:,0]))

for i in range(len(r1)):
    r1[i] = np.sqrt(np.dot(x1[i],x1[i]) + np.dot(y1[i],y1[i]) + np.dot(z1[i],z1[i]))
    r2[i] = np.sqrt(np.dot(x2[i],x2[i]) + np.dot(y2[i],y2[i]) + np.dot(z2[i],z2[i]))
    
    arg = (x1[i]*x2[i] + y1[i]*y2[i] + z1[i]*z2[i])
    theta[i] = np.arccos(arg/(r1[i]*r2[i]))
    uniform[i] = np.cos(theta[i])


#%% Bins for histogram
step_bin = 0.1
start_bin = -1
stopp_bin = 1 + step_bin
b = np.arange(start_bin, stopp_bin, step_bin)

fig, ax = plt.subplots()
ax.hist(uniform, bins=b, density=True)

#%% Bins for histogram
step_bin = 0.1
start_bin = 0
stopp_bin = np.pi + step_bin
b = np.arange(start_bin, stopp_bin, step_bin)

fig, ax = plt.subplots()
ax.hist(theta, bins=b, density=True)

#%% Centeral mean field approx. from Hartree
Z_unscreened = 2
Z_optimized  = 27/16

r = np.linspace(0, 5, 1000)
def rho(r, Z):
    return Z**3*4*r**2*np.exp(-2*Z*r)

## Bins for histogram
step_bin = 0.05
start_bin = 0
stopp_bin = 5 + step_bin
b = np.arange(start_bin, stopp_bin, step_bin)

# NORMALIZE!?!
rho_unscreened = rho(r, Z_unscreened)
rho_optimized = rho(r, Z_optimized)
fig, ax = plt.subplots()
ax.plot(r, rho_unscreened, '--', color='blue', label=r'$\rho_{unscreened}$')
ax.plot(r, rho_optimized , '--', color='red', label=r'$\rho_{screened}$')
ax.hist(r2, bins=b, density=True)

ax.set_title(r'Probability of finding a electron', fontsize='16')
ax.set_xlabel(r'$r$ [$Å$]', fontsize='16')
ax.set_ylabel(r'$\rho$', fontsize='16')

ax.grid()
ax.legend(fontsize='16', loc='upper right')

plt.savefig('../T1/figure/probability_benchmark.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%%

E = (1/N)*np.sum(El)

