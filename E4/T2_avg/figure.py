#%%
import numpy as np
import matplotlib.pyplot as plt

#%%
array = np.genfromtxt('../T2_avg/spectrum.csv', delimiter=',')
vfft = array[:,0]
freq = array[:,1]  

fig, ax = plt.subplots()
ax.plot(freq, vfft)

ax.set_xlim([0,10])



#%%
array = np.genfromtxt('../T2_avg/data/spectrum_high.csv', delimiter=',')
vfft_high = array[:,0]
freq = array[:,1]  

array = np.genfromtxt('../T2_avg/data/spectrum_low.csv', delimiter=',')
vfft_low = array[:,0]

fig, ax = plt.subplots()
ax.plot(freq, vfft_high)
ax.plot(freq, vfft_low)

ax.set_xlim([0,10])
