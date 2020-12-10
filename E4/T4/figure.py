#%%
import numpy as np
import matplotlib.pyplot as plt

#%%
array = np.genfromtxt('../T4/correlation.csv', delimiter=',')
spectrum = array[:,0]
freq = array[:,1]  

fig, ax = plt.subplots()
ax.plot(freq, spectrum)

ax.set_xlim([0,1])

#%%
array = np.genfromtxt('../T4/power_corr.csv', delimiter=',')
spectrum = array[:,0]
freq = array[:,1]

fig, ax = plt.subplots()
ax.plot(freq, spectrum)

ax.set_xlim([0,10])

