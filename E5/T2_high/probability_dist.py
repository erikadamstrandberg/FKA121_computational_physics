import numpy as np
import matplotlib.pyplot as plt

#%%
ntrails = 20

data = np.genfromtxt('../T2_high/data/timetrail0.csv', delimiter=',')

time = data[:,2]
N = len(time)
mean_x = np.zeros(N)
mean_v = np.zeros(N)
sigma2_x = np.zeros(N)
sigma2_v = np.zeros(N)

fig, ax = plt.subplots(2,1)
for i in range(ntrails):
    data = np.genfromtxt('../T2_high/data/timetrail' + str(i) + '.csv', delimiter=',')

    x = data[:,0]
    v = data[:,1]
    
    for t in range(N):
        mean_x[t] += np.sum(x[0:t])/(t+1)
        mean_v[t] += np.sum(v[0:t])/(t+1) 

    ax[0].plot(time, x, 'r', alpha=0.2)
    ax[1].plot(time, v, 'r', alpha=0.2)

mean_x = mean_x/ntrails
mean_v = mean_v/ntrails

for i in range(ntrails):
    data = np.genfromtxt('../T2_high/data/timetrail' + str(i) + '.csv', delimiter=',')

    x = data[:,0]
    v = data[:,1]
    
    for t in range(N):
        sigma2_x[t] += np.sum((x[0:t]-mean_x[0:t])**2)/(t+1)
        sigma2_v[t] += np.sum((v[0:t]-mean_v[0:t])**2)/(t+1)

sigma2_x = sigma2_x/ntrails
sigma2_v = sigma2_v/ntrails


ax[0].plot(time, mean_x, 'k')
ax[0].plot(time, mean_x + np.sqrt(sigma2_x), '--k')
ax[0].plot(time, mean_x - np.sqrt(sigma2_x), '--k')

ax[1].plot(time, mean_v, 'k')
ax[1].plot(time, mean_v + np.sqrt(sigma2_v), '--k')
ax[1].plot(time, mean_v - np.sqrt(sigma2_v), '--k')

plt.show()
