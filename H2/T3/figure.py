#%% T! Benchmarking for variational Monte Carlo
# 
# Standard includes
import numpy as np
import matplotlib.pyplot as plt

#%% Uncorrelated! With only one electrion wavefunctions

array = np.genfromtxt('../T2/local_energy.csv', delimiter=',')

El = array
fig, ax = plt.subplots()
ax.plot(El)

#%%
array = np.genfromtxt('../T2/correlation.csv', delimiter=',')

fig, ax = plt.subplots()
ax.plot(array)

#%%
array = np.genfromtxt('../T2/block_average.csv', delimiter=',')

fig, ax = plt.subplots()
ax.plot(array)