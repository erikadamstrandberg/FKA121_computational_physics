#%%
import numpy as np
import matplotlib.pyplot as plt

#%%

array = np.genfromtxt('../T2/timetrail.csv', delimiter=',')

start = 0
stop = 100
x = array[:,0]
v = array[:,1]
t = array[:,2]
N = len(t)
dt = t[1]-t[0]

fig, ax = plt.subplots(1,2)
ax[0].plot(t, x, label="x")
ax[1].plot(t, v, label="v")
plt.show()
      
array = np.genfromtxt('../T2/spectrum.csv', delimiter=',')
vfft = array[:,0]
freq = array[:,1]  

fig, ax = plt.subplots()
ax.plot(freq, vfft)

ax.set_xlim([0,10])




#%%
array = np.genfromtxt('../T2/sample.csv', delimiter=',')

fig, ax = plt.subplots()
ax.plot(array, label="x")

