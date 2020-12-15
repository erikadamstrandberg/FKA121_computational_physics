import numpy as np
import matplotlib.pyplot as plt

#%%
ntrails = 10

data = np.genfromtxt('../T2_low/data/timetrail0.csv', delimiter=',')

time = data[:,2]
N = len(time)
mean_x = np.zeros(N)
mean_v = np.zeros(N)
sigma2_x = np.zeros(N)
sigma2_v = np.zeros(N)

fig, ax = plt.subplots(2,1)
for i in range(ntrails):
    data = np.genfromtxt('../T2_low/data/timetrail' + str(i) + '.csv', delimiter=',')

    x = data[:,0]
    v = data[:,1]
    
    mean_x += x
    mean_v += v
#    for t in range(N):
#       mean_x[t] += np.sum(x[0:t])/(t+1)
#        mean_v[t] += np.sum(v[0:t])/(t+1) 

    ax[0].plot(time, x, 'r', alpha=0.2)
    ax[1].plot(time, v, 'r', alpha=0.2)

#mean_x = mean_x/ntrails
#mean_v = mean_v/ntrails
#
#for i in range(ntrails):
#    data = np.genfromtxt('../T2_low/data/timetrail' + str(i) + '.csv', delimiter=',')
#
#    x = data[:,0]
#    v = data[:,1]
#    
##    for t in range(N):
##        sigma2_x[t] += np.sum((x[0:t]-mean_x[0:t])**2)/(t+1)
##        sigma2_v[t] += np.sum((v[0:t]-mean_v[0:t])**2)/(t+1)
#    sigma2_x += (x-mean_x)**2
#    sigma2_v += (v-mean_v)**2

#sigma2_x = sigma2_x/ntrails
#sigma2_v = sigma2_v/ntrails
#print(f'Mean x: {mean_x[-1]} um. std_x: {np.sqrt(sigma2_x[-1])} um')
#print(f'Mean v: {mean_v[-1]} mm/s. std_v: {np.sqrt(sigma2_v[-1])} mm/s')

var_mean = np.genfromtxt('../T2_low/data/mean_and_var.csv', delimiter=',')

mean_x = var_mean[:,0]
mean_v = var_mean[:,1]
sigma2_x = var_mean[:,2]
sigma2_v = var_mean[:,3]

ax[0].plot(time, mean_x, 'k')
ax[0].plot(time, mean_x + np.sqrt(sigma2_x), '--k')
ax[0].plot(time, mean_x - np.sqrt(sigma2_x), '--k')
ax[0].set_title('x')
ax[0].set_xlabel('t [ms]')
ax[0].set_ylabel('x [um]')

ax[1].plot(time, mean_v, 'k')
ax[1].plot(time, mean_v + np.sqrt(sigma2_v), '--k')
ax[1].plot(time, mean_v - np.sqrt(sigma2_v), '--k')
ax[1].set_title('v')
ax[1].set_xlabel('t [ms]')
ax[1].set_ylabel('x [mm/s]')


plt.show()
