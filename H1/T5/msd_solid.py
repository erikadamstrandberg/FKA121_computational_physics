import numpy as np
import matplotlib.pyplot as plt

#%%
pos_data = np.genfromtxt('../T3/data/merged_q_trail.csv', delimiter=',')

M, n_particles = np.shape(pos_data)
time = np.arange(0,M-1,1) #pos_data[0:-1,-1]
n_particles = int(n_particles/3)
print("Number of timesteps: " + str(M))
print("Number of particles: " + str(n_particles))

msd_x = np.zeros(len(time))
msd_y = np.zeros(len(time))
msd_z = np.zeros(len(time))

for n in range(n_particles):
    x = pos_data[:,n*3]
    y = pos_data[:,n*3+1]
    z = pos_data[:,n*3+2]
    for k in range(len(time)): 
        x_forward = x[k:M]
        y_forward = y[k:M]
        z_forward = z[k:M]
        x_lag = x[0:M-k]
        y_lag = y[0:M-k]
        z_lag = z[0:M-k]
        msd_x[k] += np.sum((x_forward-x_lag)**2)/len(x_forward)
        msd_y[k] += np.sum((y_forward-y_lag)**2)/len(y_forward)
        msd_z[k] += np.sum((z_forward-z_lag)**2)/len(z_forward)


msd_x = msd_x/n_particles
msd_y = msd_y/n_particles
msd_z = msd_z/n_particles

msd_tot = msd_x + msd_y + msd_z
fig, ax = plt.subplots()
ax.plot(time, msd_x, label=r"$msd_x$")
ax.plot(time, msd_y, label=r"$msd_y$")
ax.plot(time, msd_z, label=r"$msd_z$")
ax.plot(time, msd_z+msd_y+msd_x, label=r"$msd_r$")
ax.legend()
plt.show()
#plt.savefig('figure/temp_equil_T500.pdf', format='pdf', bbox_inches='tight')


