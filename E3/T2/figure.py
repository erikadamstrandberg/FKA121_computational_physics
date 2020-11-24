#%%
import numpy as np
import matplotlib.pyplot as plt

#%%

array = np.genfromtxt('dist_10.csv', delimiter=',', skip_header=1)

uniform = array[:,0]
sinx_dist = array[:,1]

step_bin = 0.05
start_bin = 0
stopp_bin = 1 + step_bin
b = np.arange(start_bin, stopp_bin, step_bin)

fig, ax = plt.subplots()
plt.hist(sinx_dist, bins=b) 

plt.savefig('sinx_N_10.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%%

array = np.genfromtxt('dist_100.csv', delimiter=',', skip_header=1)

uniform = array[:,0]
sinx_dist = array[:,1]

step_bin = 0.05
start_bin = 0
stopp_bin = 1 + step_bin
b = np.arange(start_bin, stopp_bin, step_bin)

fig, ax = plt.subplots()
plt.hist(sinx_dist, bins=b) 

plt.savefig('sinx_N_100.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%%

array = np.genfromtxt('dist_1000.csv', delimiter=',', skip_header=1)

uniform = array[:,0]
sinx_dist = array[:,1]

step_bin = 0.05
start_bin = 0
stopp_bin = 1 + step_bin
b = np.arange(start_bin, stopp_bin, step_bin)

fig, ax = plt.subplots()
plt.hist(sinx_dist, bins=b) 

plt.savefig('sinx_N_1000.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%%

array = np.genfromtxt('dist_10000.csv', delimiter=',', skip_header=1)

uniform = array[:,0]
sinx_dist = array[:,1]

step_bin = 0.02
start_bin = 0
stopp_bin = 1 + step_bin
b = np.arange(start_bin, stopp_bin, step_bin)

fig, ax = plt.subplots()
plt.hist(sinx_dist, bins=b) 

plt.savefig('sinx_N_10000.pdf', format='pdf', bbox_inches='tight')
plt.show()


#%%

array = np.genfromtxt('dist_100000.csv', delimiter=',', skip_header=1)

uniform = array[:,0]
sinx_dist = array[:,1]

step_bin = 0.01
start_bin = 0
stopp_bin = 1 + step_bin
b = np.arange(start_bin, stopp_bin, step_bin)

fig, ax = plt.subplots()
plt.hist(sinx_dist, bins=b) 

plt.savefig('sinx_N_10000.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%%

x = np.linspace(0.01,0.99,1000)
f = x*(1-x)
p = (np.pi/2)*np.sin(np.pi*x)
g = (2/np.pi)*x*(1-x)/np.sin(np.pi*x)

fig, ax = plt.subplots()
ax.plot(x,f)
ax.plot(x,p)
ax.plot(x,g)




