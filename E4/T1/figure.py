#%%
import numpy as np
import matplotlib.pyplot as plt


case_low = np.genfromtxt('../T1/data/timetrail_low_005.csv', delimiter=',')
case_high = np.genfromtxt('../T1/data/timetrail_high_005.csv', delimiter=',')

x_low = case_low[:,0]
v_low = case_low[:,1]
x_high = case_high[:,0]
v_high = case_high[:,1]

t = case_low[:,2]

fig, ax = plt.subplots(2,2)
ax[0,0].plot(t, x_high*1e3)
ax[1,0].plot(t, v_high)
ax[0,0].set_title("P = 99.8 kPa")

ax[0,1].plot(t, x_low*1e3)
ax[1,1].plot(t, v_low)
ax[0,1].set_title("P = 2.75 kPa")


ax[0,0].set_ylabel("position [nm]")
ax[1,0].set_ylabel("velocity [mm/s]")
ax[1,0].set_xlabel("time [ms]")
ax[1,1].set_xlabel("time [ms]")

#plt.savefig("timetrail_dt_005.pdf")
plt.show()


#%%

case_low = np.genfromtxt('../T1/timetrail.csv', delimiter=',')

fig, ax = plt.subplots()
ax.plot(case_low[:,2], case_low[:,0])