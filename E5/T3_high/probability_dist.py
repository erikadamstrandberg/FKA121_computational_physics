import numpy as np
import matplotlib.pyplot as plt

#%%
ntrails = 900

data = np.genfromtxt('../T3_high/data/timetrail0.csv', delimiter=',')

time = data[:,2]
dt = time[1]-time[0]
N = len(time)
mean_x = np.zeros(N)
mean_v = np.zeros(N)
sigma2_x = np.zeros(N)
sigma2_v = np.zeros(N)

for i in range(ntrails):
    data = np.genfromtxt('../T3_high/data/timetrail' + str(i) + '.csv', delimiter=',')

    x = data[:,0]
    v = data[:,1]
    
    for t in range(N):
        mean_x[t] += np.sum(x[0:t])/(t+1)
        mean_v[t] += np.sum(v[0:t])/(t+1) 

mean_x = mean_x/ntrails
mean_v = mean_v/ntrails

for i in range(ntrails):
    data = np.genfromtxt('../T3_high/data/timetrail' + str(i) + '.csv', delimiter=',')

    x = data[:,0]
    v = data[:,1]
    
    for t in range(N):
        sigma2_x[t] += np.sum((x[0:t]-mean_x[0:t])**2)/(t+1)
        sigma2_v[t] += np.sum((v[0:t]-mean_v[0:t])**2)/(t+1)

sigma2_x = sigma2_x/ntrails
sigma2_v = sigma2_v/ntrails

# T3 
step_bin = 0.005
start_bin = -0.2
stop_bin = 0.2 + step_bin
x_bins = np.arange(start_bin, stop_bin, step_bin)

step_bin = 0.05
start_bin = -2
stop_bin = 2 + step_bin
v_bins = np.arange(start_bin, stop_bin, step_bin)

time_steps = np.array([50, 100, 150, 300, 500])
x_hist = np.zeros([ntrails, len(time_steps)])
v_hist = np.zeros([ntrails, len(time_steps)])

for i in range(ntrails):
    for j, ts in enumerate(time_steps):
        data = np.genfromtxt('../T3_high/data/timetrail' + str(i) + '.csv', delimiter=',')

        x_hist[i,j] = data[ts,0]
        v_hist[i,j] = data[ts,1]

def gauass(x, mean, var):
    return np.exp(-(x-mean)**2/(2*var))/(np.sqrt(2*np.pi*var))

fig, ax = plt.subplots(2,1)
for i, ts in enumerate(time_steps):
    x_mean = np.sum(x_hist[:,i])/ntrails
    x_var = np.sum((x_hist[:,i]-x_mean)**2)/ntrails
    ax[0].hist(x_hist[:,i], bins=x_bins, density=True, label='t = '+str(dt*ts), alpha=0.5)
    ax[0].plot(x_bins, gauass(x_bins, x_mean, x_var), '--k')
    
    v_mean = np.sum(v_hist[:,i])/ntrails
    v_var = np.sum((v_hist[:,i]-v_mean)**2)/ntrails
    ax[1].hist(v_hist[:,i], bins=v_bins, density=True, label='t = '+str(dt*ts), alpha=0.5)
    ax[1].plot(v_bins, gauass(v_bins, v_mean, v_var), '--k')

ax[0].legend(loc='upper right')
ax[0].set_title('x')
ax[0].set_xlabel('x [um]')
ax[1].legend(loc='upper right')
ax[1].set_title('v')
ax[1].set_xlabel('v [mm/s]')
plt.savefig('../T3_high/figure/dist_high_unit.pdf', format='pdf', bbox_inches='tight')
plt.show()

