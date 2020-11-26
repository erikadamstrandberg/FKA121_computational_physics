#%%
import numpy as np
import matplotlib.pyplot as plt

#%%

array = np.genfromtxt('walkers.csv', delimiter=',', skip_header=1)

x = array[:,0]
y = array[:,1]
z = array[:,2]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x,y,z)

#%%

fig, ax = plt.subplots()
ax.plot(x)

#%%

step_bin = 0.05
start_bin = 0
stopp_bin = 1 + step_bin
b = np.arange(start_bin, stopp_bin, step_bin)

fig, ax = plt.subplots()
plt.hist(sinx_dist, bins=b) 