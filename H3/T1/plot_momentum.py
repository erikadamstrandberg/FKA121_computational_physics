#%% H3b time-independent scroooodinger

# Standard includes
import numpy as np
import matplotlib.pyplot as plt

m_prim_u  = 0.009647464077950547
hbar_prim = 0.6582845299625468 

#%%
def wave_package_position(x, x0, p0, d, hbar_prim):    
    prefactor = np.sqrt(np.sqrt(np.pi)*d)
    envelope  = np.exp(-(x - x0)**2/(2*d**2))
    plane_wave = np.exp(1j*p0*(x - x0)/hbar_prim)
    return (1/prefactor)*envelope*plane_wave

def std_from_FWHM(n, dx):
    max_value = np.max(n)
    width_limit = max_value/2
    return np.sum(n > width_limit)*dx/(2*np.sqrt(2*np.log(2)))

def p_width_anal(d, hbar_prim):
    return hbar_prim/(np.sqrt(2)*d)

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

def wave_package_momentum(p, p0, d, hbar_asup):
    prefactor = np.sqrt(d/(np.sqrt(np.pi)*hbar_asup))
    envelope  = np.exp(-d**2*(p0 - p)**2/(2*hbar_asup**2))
    return prefactor*envelope

dp        = 2*np.pi*hbar_prim/(Nx*dx)              # Delta momentum
p_nyquist = 2*np.pi*hbar_prim/(2*dx)               # Nyquivst momentum
p         = np.arange(-p_nyquist, p_nyquist, dp)   # Generating momentum "freq"

phi_p = np.fft.fftshift(np.fft.fft(phi_x))         # Shifting values from transform
n_p   = np.abs(phi_p)**2*dx**2/(2*np.pi*hbar_prim) # Probability density

print(np.sum(n_p)*dp)

phi_p_anal = wave_package_momentum(p, p0, d, hbar_prim)
n_p_anal = np.abs(phi_p_anal)**2

print(np.sum(n_p_anal)*dp)

print(std_from_FWHM(n_p, dp))
print(p_width_anal(d, hbar_prim))

k = p/hbar_prim
k0 = p0/hbar_prim

plot_p0_x = np.array([k0, k0])
y_min = -0.01
y_max = np.max(n_p)+0.1
plot_p0_y = np.array([y_min, y_max])

fig, ax = plt.subplots()
ax.plot(k, n_p, color='lightgreen', linewidth=4      , label=r'$|\psi(k)|^2$')
ax.plot(k, n_p_anal, '--', color='black', linewidth=4, label=r'$|\psi(k)|^2_{analytical}$')
ax.plot(plot_p0_x, plot_p0_y, color='black'          , label=r'$k_0=6.91\ Å^{-1}$')

ax.set_title(r'Momentum of initial Gaussian wave packet', fontsize='16')
ax.set_xlabel(r'$k$ [$Å^{-1}$]', fontsize='16')
ax.set_ylabel(r'$|\psi(k)|^2$ [$Å$]', fontsize='16')

ax.grid()
ax.legend(fontsize='16', loc='upper right')

ax.set_xlim([0,14])
ax.set_ylim([y_min,y_max])

#plt.savefig('../T1/initial_gauss_momentum.pdf', format='pdf', bbox_inches='tight')
#plt.show()
