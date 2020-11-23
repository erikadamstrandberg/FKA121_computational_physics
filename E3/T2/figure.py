#%%
import numpy as np
import matplotlib.pyplot as plt

#%%

array = np.genfromtxt('dist_10.csv', delimiter=',', skip_header=1)

uniform = array[:,0]
sinx_dist = array[:,1]

b = np.arange(0,1.1,0.05)

fig, ax = plt.subplots()
plt.hist(sinx_dist, bins=b) 

plt.savefig('sinx_N_10.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%%

array = np.genfromtxt('dist_100.csv', delimiter=',', skip_header=1)

uniform = array[:,0]
sinx_dist = array[:,1]

b = np.arange(0,1.1,0.05)

fig, ax = plt.subplots()
plt.hist(sinx_dist, bins=b) 

plt.savefig('sinx_N_100.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%%

array = np.genfromtxt('dist_1000.csv', delimiter=',', skip_header=1)

uniform = array[:,0]
sinx_dist = array[:,1]

b = np.arange(0,1.1,0.02)

fig, ax = plt.subplots()
plt.hist(sinx_dist, bins=b) 

plt.savefig('sinx_N_1000.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%%

array = np.genfromtxt('dist_10000.csv', delimiter=',', skip_header=1)

uniform = array[:,0]
sinx_dist = array[:,1]

b = np.arange(0,1.1,0.02)

fig, ax = plt.subplots()
plt.hist(sinx_dist, bins=b) 

plt.savefig('sinx_N_10000.pdf', format='pdf', bbox_inches='tight')
plt.show()


#%%

array = np.genfromtxt('dist_100000.csv', delimiter=',', skip_header=1)

uniform = array[:,0]
sinx_dist = array[:,1]

b = np.arange(0,1.1,0.001)

fig, ax = plt.subplots()
plt.hist(sinx_dist, bins=b) 

plt.savefig('sinx_N_10000.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%%


