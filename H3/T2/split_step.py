import numpy as np
import matplotlib.pyplot as plt

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

hbar_si = 1.054571817e-34           # hbar in SI
hbar_prim = hbar_si*1e34/1.602      # hbar in ASU_prim. eV*fs
m_prim_u = 0.009647                 # mass/u in ASU_prim


#%%
dt = 0.05

d   = 0.5               # Width of our hydrogen atom
m_h = 1*m_prim_u        # Mass of our hydrogen atom

p0 = np.sqrt(0.2*m_h)   # Initial momentum of our hydrogen atom
x0 = 0.0                # Initial position of our hydrogen atom

dx      = 0.001
x_start = -50.0
x_stop  = 50.0
x = np.arange(x_start, x_stop, dx)
Nx = len(x)

phi_x = wave_package_position(x, x0, p0, d, hbar_prim)
n_x = np.abs(phi_x)**2

# fourier transform

dp        = 2*np.pi*hbar_prim/(Nx*dx)              # Delta momentum
p_nyquist = 2*np.pi*hbar_prim/(2*dx)               # Nyquivst momentum
p         = np.arange(-p_nyquist, p_nyquist, dp)   # Generating momentum "freq"
p         = np.fft.fftshift(p)

v_prop = 1
p_prop = np.exp(-1j*p**2*dt/(hbar_prim*2*m_h))

phi_x_2 = propagate(phi_x, p_prop, v_prop)

phi_x_anal = wave_package_position_evo(x, x0, dt, p0, d, m_h, hbar_prim)

fig, ax = plt.subplots()
#ax.plot(x, phi_x)
ax.plot(x, np.imag(phi_x_2))
ax.plot(x, np.imag(phi_x_anal), '--')

plt.show()


