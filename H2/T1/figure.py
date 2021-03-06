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
step_bin = 0.05
start_bin = -1
stopp_bin = 1 + step_bin
b = np.arange(start_bin, stopp_bin, step_bin)

x = np.array([-1,1])
y = np.array([1/2,1/2])

fig, ax = plt.subplots()
ax.hist(uniform, bins=b, density=True, label=r'$x=\cos(\theta)$')
ax.plot(x,y, '--', color='black', label=r'$1/2$', linewidth=2)

ax.set_title(r'Distribution of $x=\cos(\theta)$ for $\Psi_{un}$', fontsize='16')
ax.set_xlabel(r'$x$', fontsize='16')
ax.set_ylabel(r'$\rho$', fontsize='16')

ax.set_ylim([0,0.8])

ax.grid()
ax.legend(fontsize='16', loc='upper right')

plt.savefig('../T1/figure/uniform_uncorrelated.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%% Bins for histogram
step_bin = 0.05
start_bin = 0
stopp_bin = np.pi + step_bin
b = np.arange(start_bin, stopp_bin, step_bin)

x = np.linspace(0,np.pi,500)
y = (1/2)*np.sin(x)


fig, ax = plt.subplots()
ax.hist(theta, bins=b, density=True, label=r'$\theta$')
ax.plot(x,y, '--', color='black', label=r'$1/2\sin(\theta)$', linewidth=2)

ax.set_title(r'Distribution of $\theta$ for $\Psi_{un}$', fontsize='16')
ax.set_xlabel(r'$\theta$', fontsize='16')
ax.set_ylabel(r'$\rho$', fontsize='16')

ax.grid()
ax.legend(fontsize='16', loc='upper left')

plt.savefig('../T1/figure/theta_uncorrelated.pdf', format='pdf', bbox_inches='tight')
plt.show()


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
step_bin = 0.05
start_bin = -1
stopp_bin = 1 + step_bin
b = np.arange(start_bin, stopp_bin, step_bin)

x = np.array([-1,1])
y = np.array([1/2,1/2])

fig, ax = plt.subplots()
ax.hist(uniform, bins=b, density=True, label=r'$x=\cos(\theta)$')
ax.plot(x,y, '--', color='black', label=r'$1/2$', linewidth=2)

ax.set_title(r'Distribution of $x=\cos(\theta)$', fontsize='16')
ax.set_xlabel(r'$x$', fontsize='16')
ax.set_ylabel(r'$\rho$', fontsize='16')

ax.grid()
ax.legend(fontsize='16', loc='upper right')

plt.savefig('../T1/figure/uniform_correlated.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%% Bins for histogram
step_bin = 0.05
start_bin = 0
stopp_bin = np.pi + step_bin
b = np.arange(start_bin, stopp_bin, step_bin)

x = np.linspace(0,np.pi,500)
y = (1/2)*np.sin(x)


fig, ax = plt.subplots()
ax.hist(theta, bins=b, density=True, label=r'$\theta$')
ax.plot(x,y, '--', color='black', label=r'$1/2\sin(\theta)$', linewidth=2)

ax.set_title(r'Distribution of $\theta$', fontsize='16')
ax.set_xlabel(r'$\theta$', fontsize='16')
ax.set_ylabel(r'$\rho$', fontsize='16')

ax.grid()
ax.legend(fontsize='16', loc='upper left')

plt.savefig('../T1/figure/theta_correlated.pdf', format='pdf', bbox_inches='tight')
plt.show()


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
ax.hist(r1, bins=b, density=True, label=r'$\rho(r_1)$')
ax.plot(r, rho_unscreened, '--', color='black', label=r'$\rho_{unscreened}(r)$', linewidth=2)
ax.plot(r, rho_optimized , '--', color='red', label=r'$\rho_{screened}(r)$', linewidth=2)

ax.set_title(r'Benchmarking $\rho(r_1)$ sampling, $\alpha=0.1$', fontsize='16')
ax.set_xlabel(r'$r$ [$a_0$]', fontsize='16')
ax.set_ylabel(r'$\rho$', fontsize='16')

ax.grid()
ax.legend(fontsize='16', loc='upper right')

plt.savefig('../T1/figure/probability_benchmark_r1.pdf', format='pdf', bbox_inches='tight')
plt.show()

fig, ax = plt.subplots()
ax.hist(r2, bins=b, density=True, label=r'$\rho(r_2)$')
ax.plot(r, rho_unscreened, '--', color='black', label=r'$\rho_{unscreened}(r)$', linewidth=2)
ax.plot(r, rho_optimized , '--', color='red', label=r'$\rho_{screened}(r)$', linewidth=2)

ax.set_title(r'Benchmarking $\rho(r_2)$ sampling, $\alpha=0.1$', fontsize='16')
ax.set_xlabel(r'$r$ [$a_0$]', fontsize='16')
ax.set_ylabel(r'$\rho$', fontsize='16')

ax.grid()
ax.legend(fontsize='16', loc='upper right')

plt.savefig('../T1/figure/probability_benchmark_r2.pdf', format='pdf', bbox_inches='tight')
plt.show()

r_tot = np.append(r1,r2)
fig, ax = plt.subplots()
ax.hist(r_tot, bins=b, density=True, label=r'$\rho(r)$')
ax.plot(r, rho_unscreened, '--', color='black', label=r'$\rho_{unscreened}(r)$', linewidth=2)
ax.plot(r, rho_optimized , '--', color='red', label=r'$\rho_{screened}(r)$', linewidth=2)

ax.set_title(r'Benchmarking $\rho(r)$ sampling, $\alpha=0.1$', fontsize='16')
ax.set_xlabel(r'$r$ [$a_0$]', fontsize='16')
ax.set_ylabel(r'$\rho$', fontsize='16')

ax.grid()
ax.legend(fontsize='16', loc='upper right')

plt.savefig('../T1/figure/probability_benchmark_r.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%%

r1 = np.arange(-1,1)
r2 = np.arange(-1,1)
phi = np.exp(-2*r1)*np.exp(-2*r2)

fig = plt.figure()
ax = fig.gca(projection='3d')

# Make data.
r1 = np.arange(0, 1, 0.01)
r2 = np.arange(0, 1, 0.01)
R1, R2 = np.meshgrid(r1, r2)
theta = 60*np.pi/180
R12 = np.sqrt(R1**2 + R2**2 - 2*R1*R2*np.cos(theta))
Z = np.exp(-2*R1)*np.exp(-2*R2)*np.exp(R12/(2*(1+0.9*R12)))
surf = ax.plot_surface(R1, R2, Z)

plt.show()

#%%

E = (1/N)*np.sum(El)



