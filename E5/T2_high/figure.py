#%%
import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('../T2/timetrail.csv', delimiter=',')

x = data[:,0]
v = data[:,1]
t = data[:,2]

fig, ax = plt.subplots()
ax.plot(t, v)

#plt.savefig("timetrail_dt_005.pdf")
plt.show()
