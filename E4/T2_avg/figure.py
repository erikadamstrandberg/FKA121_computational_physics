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
