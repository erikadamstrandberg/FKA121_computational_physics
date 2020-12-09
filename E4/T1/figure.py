#%%
import numpy as np
import matplotlib.pyplot as plt


# array[:,0] = x
# array[:,1] = v
# array[:,2] = t

array = np.genfromtxt('../T1/timetrail.csv', delimiter=',')

start = 0
stop = 100
x = array[:,0]
v = array[:,1]
t = array[:,2]

fig, ax = plt.subplots(1,2)
ax[0].plot(t, x, label="x")
ax[1].plot(t, v, label="v")
plt.show()
