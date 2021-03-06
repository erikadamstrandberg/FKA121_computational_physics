import numpy as np
import matplotlib.pyplot as plt

#%%
# STOKES LAW Fd = 6 pi mu R v
# from langevin eta = 6 pi mu R

mu = 18.6e-6
R = 2.79e-6/2
V = 4*np.pi*R**3/3
rho = 2650
m = rho*V

tau_high = 147.3e-6
tau_low  = 48.5e-6
mu_high = m/tau_high
mu_low = m/tau_low

stokes_term = 6*np.pi*mu*R

print(f'Stokes term: {stokes_term}')
print(f'Friction terms in langvin. Case high: {mu_high}. Case low: {mu_low}')


#%%
ntrails = 99

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
    
#    mean_x += x
#    mean_v += v
#
    ax[0].plot(time, x, 'r', alpha=0.2)
    ax[1].plot(time, v, 'r', alpha=0.2)
#
#mean_x = mean_x/ntrails
#mean_v = mean_v/ntrails

#for i in range(ntrails):
#    data = np.genfromtxt('../T2_high/data/timetrail' + str(i) + '.csv', delimiter=',')
#
#    x = data[:,0]
#    v = data[:,1]
#    
#    sigma2_x += (x-mean_x)**2
#    sigma2_v += (v-mean_v)**2
#
#sigma2_x = sigma2_x/ntrails
#sigma2_v = sigma2_v/ntrails
    
var_mean = np.genfromtxt('../T2_high/data/mean_and_var.csv', delimiter=',')

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
