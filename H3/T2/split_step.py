#%% H3b time-independent scroooodinger

import numpy as np
import matplotlib.pyplot as plt

hbar_si = 1.054571817e-34           # hbar in SI
hbar_prim = hbar_si*1e34/1.602      # hbar in ASU_prim. eV*fs
m_prim_u = 0.009647464077950547     # mass/u in ASU_prim

#%%

def wave_package_position_evo(x, x0, t, p0, d, m, hbar_prim):   
    prefactor = np.sqrt(np.sqrt(np.pi)*(d+1j*hbar_prim*t/(m*d)))
    envelope  = np.exp(-(x - x0 - p0*t/m)**2/(2*d**2*(1 + 1j*hbar_prim*t/(m*d**2))))
    plane_wave = np.exp(1j*p0*(x - x0 - p0*t/(2*m))/hbar_prim)
    return (1/prefactor)*envelope*plane_wave

def wave_package_position(x, x0, p0, d, hbar_prim):    
    prefactor = np.sqrt(np.sqrt(np.pi)*d)
    envelope  = np.exp(-(x - x0)**2/(2*d**2))
    plane_wave = np.exp(1j*p0*(x - x0)/hbar_prim)
    return (1/prefactor)*envelope*plane_wave

def propagate(phi_x, p_prop, v_prop):
    return np.fft.ifft(p_prop*np.fft.fft(v_prop*phi_x))

def std_from_FWHM(n, dx):
    max_value = np.max(n)
    width_limit = max_value/2
    return np.sum(n > width_limit)*dx/(2*np.sqrt(2*np.log(2)))

def x_width_anal(t, d, m, hbar_prim):
    prefactor = d/np.sqrt(2)
    time_factor = np.sqrt(1 + hbar_prim**2*t**2/(m**2*d**4))
    return prefactor*time_factor

def p_width_anal(t, d, hbar_prim, m):
    return hbar_prim/(np.sqrt(2)*d)

#%% Generate initial wave packet

dt = 0.05
d   = 0.5               # Width of our hydrogen atom
m_h = 1/m_prim_u        # Mass of our hydrogen atom

p0 = np.sqrt(0.2*m_h)   # Initial momentum of our hydrogen atom
x0 = 0.0                # Initial position of our hydrogen atom

dx      = 0.001
x_start = -50.0
x_stop  = 50.0
x = np.arange(x_start, x_stop, dx)
Nx = len(x)

phi_x = wave_package_position(x, x0, p0, d, hbar_prim)

## Set up split-step propagation
T = 420
dt = 0.5
Nt = int(T/dt)
time = np.arange(0,T,dt)


## Momentum
dp        = 2*np.pi*hbar_prim/(Nx*dx)              # Delta momentum
p_nyquist = 2*np.pi*hbar_prim/(2*dx)               # Nyquivst momentum
p         = np.arange(-p_nyquist, p_nyquist, dp)   # Generating momentum "freq"
p         = np.fft.fftshift(p)

## Propagators for the split-step
v_prop = 1
p_prop = np.exp(-1j*p**2*dt/(hbar_prim*2*m_h))


## Saving the probabilty densities and widths
save_every = np.array([100, 200, 300, 400])
n_x_plot = np.zeros([len(save_every)+1, Nx])
n_p_plot = np.zeros([len(save_every)+1, Nx])
n_x_plot[0] = np.abs(phi_x)**2
n_p_plot[0] = np.abs(np.fft.fftshift(np.fft.fft(phi_x)))**2

width_x = np.zeros(Nt)
width_x_anal = np.zeros(Nt)
width_p = np.zeros(Nt)
width_p_anal = np.zeros(Nt)

count = 1
for t in range(Nt):
    ## Calculate the probability densities for position and momentum
    n_x = np.abs(phi_x)**2
    phi_p = np.fft.fftshift(np.fft.fft(phi_x))
    n_p = np.abs(phi_p)**2
    
    ## Finding the position width from FWHM and analytcally
    width_x[t] = std_from_FWHM(n_x, dx)
    width_x_anal[t] = x_width_anal(time[t], d, m_h, hbar_prim)
    
    ## Finding the position width from FWHM and analytcally
    width_p[t] = std_from_FWHM(n_p, dp)
    width_p_anal[t] = p_width_anal(t, d, hbar_prim, m_h)#p_width_anal(d, hbar_prim)
    
    ## Save n_x for plotting
    if t*dt  in save_every:
        n_x_plot[count] = n_x
        n_p_plot[count] = n_p
        count += 1
        
    ## Propagation with the split-step
    phi_x = propagate(phi_x, p_prop, v_prop)
    

phi_x_initial_anal = wave_package_position_evo(x, x0, 0, p0, d, m_h, hbar_prim)
n_x_initial_anal   = np.abs(phi_x_initial_anal)**2


#%%
fig, ax = plt.subplots()
ax.plot(time, width_x, color='lightblue', linewidth=4, label=r'$\sigma_x(t)$')
ax.plot(time, width_x_anal, '--', color='black', linewidth=4, label=r'$\sigma_{x,analytical}(t)$')

ax.plot(time, width_p/hbar_prim, color='red', linewidth=4, label=r'$\sigma_k(t)$')
ax.plot(time, width_p_anal/hbar_prim, '--', color='black', linewidth=4, label=r'$\sigma_{k,analytical}(t)$')

ax.set_title(r'Time evolution of width for Gaussian wave packet', fontsize='16')
ax.set_xlabel(r'$t$ [$fs$]', fontsize='16')
ax.set_ylabel(r'$x$ [$Å$], $k$ [$Å^{-1}$]', fontsize='16')

ax.grid()
ax.legend(fontsize='16', loc='upper left')

plt.savefig('../T2/evolution_of_width.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%%
fig, ax = plt.subplots()

for i in range(len(save_every) +1):
    if i == 0:
        ax.plot(x, n_x_plot[i], linewidth=4, label=r't = 0 fs'.format(save_every[i]))
    else:
        ax.plot(x, n_x_plot[i], linewidth=4, label=r't = {:.0f} fs'.format(save_every[i-1]))

for i in range(len(save_every)):
    phi_x_anal = wave_package_position_evo(x, x0, save_every[i], p0, d, m_h, hbar_prim)
    n_x_anal   = np.abs(phi_x_anal)**2
    ax.plot(x, n_x_anal, '--', color='black', linewidth=2)

ax.plot(x, n_x_initial_anal, '--', color='black', linewidth=4)

ax.set_title(r'Free space propagation', fontsize='16')
ax.set_xlabel(r'$x$ [$Å$]', fontsize='16')
ax.set_ylabel(r'$|\psi(k)|^2$ [$Å^{-1}$]', fontsize='16')

ax.grid()
ax.legend(fontsize='16', loc='upper right')

ax.set_xlim([-2,30])

plt.savefig('../T2/evolution_wave_packet.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%%

k = np.fft.ifftshift(p)/hbar_prim

fig, ax = plt.subplots()
for i in range(len(save_every) + 1):
    ax.plot(k, n_p_plot[i])


ax.set_xlim([-0,20])

