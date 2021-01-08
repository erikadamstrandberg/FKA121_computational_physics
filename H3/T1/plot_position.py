#%% H3b time-independent scroooodinger

# Standard includes
import numpy as np
import matplotlib.pyplot as plt

m_prim_u  = 0.009647464077950547
hbar_prim = 0.6582845299625468 

#%% Plotting the Gaussian wave packet

def wave_package_position(x, x0, p0, d, hbar_prim):    
    prefactor = np.sqrt(np.sqrt(np.pi)*d)
    envelope  = np.exp(-(x - x0)**2/(2*d**2))
    plane_wave = np.exp(1j*p0*(x - x0)/hbar_prim)
    return (1/prefactor)*envelope*plane_wave

def wave_package_position_evo(x, x0, t, p0, d, m, hbar_prim):   
    prefactor = np.sqrt(np.sqrt(np.pi)*(d+1j*hbar_prim*t/(m*d)))
    envelope  = np.exp(-(x - x0 - p0*t/m)**2/(2*d**2*(1 + 1j*hbar_prim*t/(m*d**2))))
    plane_wave = np.exp(1j*p0*(x - x0 - p0*t/(2*m))/hbar_prim)
    return (1/prefactor)*envelope*plane_wave

def x_width_anal(t, d, m, hbar_prim):
    prefactor = d/np.sqrt(2)
    time_factor = np.sqrt(1 + hbar_prim**2*t**2/(m**2*d**4))
    return prefactor*time_factor

def std_from_FWHM(n, dx):
    max_value = np.max(n)
    width_limit = max_value/2
    return np.sum(n > width_limit)*dx/(2*np.sqrt(2*np.log(2)))

d   = 0.5               # Width of our hydrogen atom
m_h = 1/m_prim_u        # Mass of our hydrogen atom

p0 = np.sqrt(0.2*m_h)  # Initial momentum of our hydrogen atom
x0 = 3.0               # Initial position of our hydrogen atom

dx      = 0.0001
x_start = -50
x_stop  = 50.0
x       = np.arange(x_start, x_stop, dx)
Nx      = len(x)

phi_x = wave_package_position(x, x0, p0, d, hbar_prim)
n_x   = np.abs(phi_x)**2

t1      = 20
phi_x_1 = wave_package_position_evo(x, x0, t1, p0, d, m_h, hbar_prim)
n_x_1   = np.abs(phi_x_1)**2

width = std_from_FWHM(n_x_1, dx)
width_anal = x_width_anal(t1, d, m_h, hbar_prim)

total_probability = np.sum(n_x_1*dx)

print(f'Width from finding the FWHM of the wavepacket: {width} Å')
print(f'Analytical width:                              {width_anal} Å')
print(f'Total probability: {total_probability}')

plot_x0_x = np.array([x0,x0])
y_min = -0.01
y_max = np.max(n_x)+0.1
plot_x0_y = np.array([y_min,y_max])

fig, ax = plt.subplots()
ax.plot(x, phi_x)


#%%


fig, ax = plt.subplots()
ax.plot(x, n_x, color='blue', linewidth=4, label=r'$|\psi(x)|^2$')
ax.plot(plot_x0_x, plot_x0_y, color='black', label=r'$x_0=3.0\ Å$')

ax.set_title(r'Position of initial Gaussian wave packet', fontsize='16')
ax.set_xlabel(r'$x$ [$Å$]', fontsize='16')
ax.set_ylabel(r'$|\psi(x)|^2$ [$Å^{-1}$]', fontsize='16')

ax.grid()
ax.legend(fontsize='16', loc='upper right')

ax.set_xlim([0,6])
ax.set_ylim([y_min,y_max])

plt.savefig('../T1/initial_gauss_position.pdf', format='pdf', bbox_inches='tight')
plt.show()