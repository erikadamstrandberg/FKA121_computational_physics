
import numpy as np
import matplotlib.pyplot as plt


#%%
vel_data = np.genfromtxt('../T3/v_trail.csv', delimiter=',') # Å/ps

M, n_particles = np.shape(vel_data)
time = vel_data[0:-1,-1]
n_particles = int(n_particles/3)
print("Number of saved timesteps: " + str(M))

# calculate phi(0)
m = 26/9649     
kb  = 8.617e-5  # eV/K

phi_0 = 0.0
for i in range(n_particles):
    phi_0 += np.sum(vel_data[:,i]**2)/M
phi_0 = phi_0/n_particles

# Calculate velocity correlation
phi_x = np.zeros(len(time))
phi_y = np.zeros(len(time))
phi_z = np.zeros(len(time))

for n in range(n_particles):
    x = vel_data[:,n*3]
    y = vel_data[:,n*3+1]
    z = vel_data[:,n*3+2]
    for k in range(len(time)): 
        x_forward = x[k:M]
        y_forward = y[k:M]
        z_forward = z[k:M]
        x_lag = x[0:M-k]
        y_lag = y[0:M-k]
        z_lag = z[0:M-k]
        phi_x[k] += np.dot(x_forward, x_lag)/len(x_forward)
        phi_y[k] += np.dot(y_forward, y_lag)/len(y_forward)
        phi_z[k] += np.dot(z_forward, z_lag)/len(z_forward)

phi_x = phi_x/n_particles
phi_y = phi_y/n_particles
phi_z = phi_z/n_particles

phi_tot = phi_x + phi_y + phi_z

fig, ax = plt.subplots()
ax.plot(time, phi_tot, label=r"$\phi(t)$")
ax.legend()
plt.show()
#plt.savefig('figure/temp_equil_T500.pdf', format='pdf', bbox_inches='tight')


#%%
vel_data = np.genfromtxt('../T4/v_trail.csv', delimiter=',') # Å/ps

M, n_particles = np.shape(vel_data)
time = vel_data[0:-1,-1]
n_particles = int(n_particles/3)
print("Number of saved timesteps: " + str(M))

# calculate phi(0)
m = 26/9649     
kb  = 8.617e-5  # eV/K

phi_0 = 0.0
for i in range(3*n_particles):
    phi_0 += np.sum(vel_data[:,i]**2)/M
phi_0 = phi_0/n_particles
    
print("phi_0 = " + str(phi_0))

# Calculate velocity correlation
phi_x = np.zeros(len(time))
phi_y = np.zeros(len(time))
phi_z = np.zeros(len(time))

for n in range(n_particles):
    x = vel_data[:,n*3]
    y = vel_data[:,n*3+1]
    z = vel_data[:,n*3+2]
    for k in range(len(time)): 
        x_forward = x[k:M]
        y_forward = y[k:M]
        z_forward = z[k:M]
        x_lag = x[0:M-k]
        y_lag = y[0:M-k]
        z_lag = z[0:M-k]
        phi_x[k] += np.dot(x_forward, x_lag)/len(x_forward)
        phi_y[k] += np.dot(y_forward, y_lag)/len(y_forward)
        phi_z[k] += np.dot(z_forward, z_lag)/len(z_forward)

phi_x = phi_x/n_particles
phi_y = phi_y/n_particles
phi_z = phi_z/n_particles

phi_tot = phi_z + phi_y + phi_x
fig, ax = plt.subplots()
ax.plot(time, phi_tot/phi_0, label=r"$\phi(t)$")
ax.legend()
plt.show()
#plt.savefig('figure/temp_equil_T500.pdf', format='pdf', bbox_inches='tight')

#%%

delta_t = time[1] - time[0]
delta_f = 1/(time[-1])
f_max = 1/(2*delta_t)

f = np.arange(0, f_max, delta_f)
omega = 2*np.pi*f

phi_spectrum = np.zeros(len(omega))
for i in range(len(omega)):
    cos = np.dot(phi_tot, np.cos(time*omega[i]))
    phi_spectrum[i] = 2*np.sum(cos)/M
    
Ds = phi_spectrum[0]/6
print(Ds)
fig, ax = plt.subplots()
ax.plot(omega, phi_spectrum, label=r"$\phi(t)$")

#%%


