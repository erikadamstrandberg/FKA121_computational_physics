#%%
import numpy as np
import matplotlib.pyplot as plt


#%%

array = np.genfromtxt('../T3/data/corr_high.csv', delimiter=',')
phi_high = array[:,0]
t_high = array[:,1]

index_1ms = np.sum(t_high < 1)
phi_high_1ms = phi_high[0:index_1ms]
t_high_1ms = t_high[0:index_1ms]

array = np.genfromtxt('../T3/data/corr_low.csv', delimiter=',')
phi_low = array[:,0]
t_low = array[:,1]

index_1ms = np.sum(t_low < 1)
phi_low_1ms = phi_low[0:index_1ms]
t_low_1ms = t_low[0:index_1ms]

fig, ax = plt.subplots()
ax.plot(t_high_1ms, phi_high_1ms/np.max(phi_high_1ms), label="x")
ax.plot(t_low_1ms, phi_low_1ms/np.max(phi_low_1ms), label="x")

