
import numpy as np
import matplotlib.pyplot as plt



vel_data = np.genfromtxt('../T4/q_trail.csv', delimiter=',', skip_header=1) # Å/ps

M, n_particles = np.shape(vel_data)
n_particles = int(n_particles/3)
print("Number of saved timesteps: " + str(M))

# calculate phi(0)
m = 26/9649     
kb  = 8.617e-5  # eV/K

phi_0 = np.sum(vel_data**2)/M
print("phi_0 = " + str(phi_0))

# Calculate velocity correlation
t_ref = np.arange(1,500,1)  # k = t/dt
phi_x = np.zeros(len(t_ref))
phi_y = np.zeros(len(t_ref))
phi_z = np.zeros(len(t_ref))

for n in range(n_particles):
    x = vel_data[:,n*3]
    y = vel_data[:,n*3+1]
    z = vel_data[:,n*3+2]
    for i, k in enumerate(t_ref): 
        x_forward = x[k:M]
        y_forward = y[k:M]
        z_forward = z[k:M]
        x_lag = x[0:M-k]
        y_lag = y[0:M-k]
        z_lag = z[0:M-k]
        phi_x[i] += np.dot(x_forward, x_lag)/len(x_forward)
        phi_y[i] += np.dot(y_forward, y_lag)/len(y_forward)
        phi_z[i] += np.dot(z_forward, z_lag)/len(z_forward)
    phi_x = phi_x/n_particles
    phi_y = phi_y/n_particles
    phi_z = phi_z/n_particles

fig, ax = plt.subplots()
ax.plot(t_ref, phi_z + phi_y + phi_x, label=r"$\phi(t)$")
ax.legend()
plt.show()
#plt.savefig('figure/temp_equil_T500.pdf', format='pdf', bbox_inches='tight')
