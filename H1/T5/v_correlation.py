import numpy as np
import matplotlib.pyplot as plt

# velocity correlation for solid phase
vel_data = np.genfromtxt('../T3/data/merged_v_trail.csv', delimiter=',') # Å/ps

M, n_particles = np.shape(vel_data)
n_particles = int(n_particles/3)

time = vel_data[:,-1]
dt = time[1]-time[0]
time = np.arange(0,len(time)*dt,dt)

print("Number of saved timesteps: " + str(M))
print("Number of particles: " + str(n_particles))

phi_0_solid = 0.0
for i in range(3*n_particles):
    phi_0_solid += np.sum(vel_data[:,i]**2)/M
phi_0_solid = phi_0_solid/n_particles


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

phi_tot_solid = (phi_x + phi_y + phi_z)/phi_0_solid



# velocity correlation for liquid phase
vel_data = np.genfromtxt('../T4/data/merged_v_trail.csv', delimiter=',') # Å/ps

M, n_particles = np.shape(vel_data)
n_particles = int(n_particles/3)

time = vel_data[:,-1]
dt = time[1]-time[0]
time = np.arange(0,len(time)*dt,dt)

print("Number of saved timesteps: " + str(M))
print("Number of particles: " + str(n_particles))

phi_0_liquid = 0.0
for i in range(3*n_particles):
    phi_0_liquid += np.sum(vel_data[:,i]**2)/M
phi_0_liquid = phi_0_liquid/n_particles


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

phi_tot_liquid = (phi_x + phi_y + phi_z)/phi_0_liquid

# Plot velocity correlation for solid and liquid phase
fig, ax = plt.subplots()

ax.plot(time, phi_tot_solid, color='blue', label=r'solid $\Phi(t)$')
ax.plot(time, phi_tot_liquid, color='red', label=r'liquid $\Phi(t)$')

ax.set_xlim([0,1])
ax.set_xlabel(r'$t$ [$ps$]', fontsize=16)
ax.set_ylabel(r'$\Phi/\Phi_0$', fontsize=16)
ax.legend(fontsize=16)
ax.set_title(r'Velocity correlation for solid and liquid phase', fontsize=16)
ax.grid()
plt.savefig('velocity_correlation.pdf', format='pdf', bbox_inches='tight')
plt.show()


# Calculte spectrum for solid and liquid phase
delta_t = time[1] - time[0]
delta_f = 1/(M*dt)
f_max = 1/(2*delta_t)

f = np.arange(0, f_max, delta_f)
omega = 2*np.pi*f

# calculate spectrum for solid
phi_spectrum_solid = np.zeros(len(omega))
for i in range(len(omega)):
    cos = np.dot(phi_tot_solid*phi_0_solid, np.cos(time*omega[i]))
    phi_spectrum_solid[i] = 2*cos/M

# calculate spectrum for liquid 
phi_spectrum_liquid = np.zeros(len(omega))
for i in range(len(omega)):
    cos = np.dot(phi_tot_liquid*phi_0_liquid, np.cos(time*omega[i]))
    phi_spectrum_liquid[i] = 2*cos/M
    
Ds = phi_spectrum_liquid[0]/6
print("Ds: "  + str(Ds))

# plot spectrum for solid and liquid phase
fig, ax = plt.subplots()
ax.plot(omega/(2*np.pi), phi_spectrum_solid, color='blue', label=r"solid $\hat{\phi}(f)$")
ax.plot(omega/(2*np.pi), phi_spectrum_liquid, color='red', label=r"liquid $\hat{\phi}(f)$")

ax.set_xlim([0,50])
ax.set_xlabel(r'$f$ [$ps^{-1}$]', fontsize=16)
ax.set_ylabel(r'$\hat{\Phi}$ [$Å^2/ps$]', fontsize=16)
ax.legend(fontsize=16)
ax.set_title(r'Spectrum of velocity correlation  for solid and liquid phase', fontsize=16)
ax.grid()
plt.savefig('vcorr_spectrum.pdf', format='pdf', bbox_inches='tight')
plt.show()


