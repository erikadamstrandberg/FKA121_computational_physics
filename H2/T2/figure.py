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
N = len(array)
N = np.arange(0,N)

limit_for_ns = 0.135

x = np.array([np.min(N), np.max(N)])
y = np.array([limit_for_ns, limit_for_ns])

what_N_index = np.sum(array > limit_for_ns)
print(f'ns = {what_N_index}')

fig, ax = plt.subplots()
ax.plot(N, array)
ax.plot(x, y)

#%%
array = np.genfromtxt('../T2/block_average.csv', delimiter=',')

fig, ax = plt.subplots()
ax.plot(array)