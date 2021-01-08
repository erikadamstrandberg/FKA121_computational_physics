import numpy as np
import matplotlib.pyplot as plt

def wave_package_position(x, x0, p0, d, hbar_prim):    
    prefactor = np.sqrt(np.sqrt(np.pi)*d)
    envelope  = np.exp(-(x - x0)**2/(2*d**2))
    plane_wave = np.exp(1j*p0*(x - x0)/hbar_prim)
    return (1/prefactor)*envelope*plane_wave

def propagate(phi_x, p_prop, v_prop):
    return np.fft.ifft(p_prop*np.fft.fft(v_prop*phi_x))

#%% define constants

hbar_si = 1.054571817e-34           # hbar in SI
hbar_prim = hbar_si*1e34/1.602      # hbar in ASU_prim. eV*fs
m_prim_u = 0.009647                 # mass/u in ASU_prim


# position space
dx      = 0.001
x_start = -30.0
x_stop  = 30.0
x = np.arange(x_start, x_stop, dx)
Nx = len(x)

# momentum space
dp        = 2*np.pi*hbar_prim/(Nx*dx)              # Delta momentum
p_nyquist = 2*np.pi*hbar_prim/(2*dx)               # Nyquivst momentum
p         = np.arange(-p_nyquist, p_nyquist, dp)   # Generating momentum "freq"
p         = np.fft.fftshift(p)

# initial values for wavefunction
d   = 0.5               # Width of our hydrogen atom
m_h = 1/m_prim_u        # Mass of our hydrogen atom

initial_energy = 0.12
p0 = np.sqrt(2*initial_energy*m_h)   # Initial momentum of our hydrogen atom
x0 = -15                 # Initial position of our hydrogen atom

#%% Propagation
T = 550
dt = 0.1
Nt = int(T/dt)

V0 = 0.1
alpha = 2.0

V_x = V0/(np.cosh(x/alpha)**2)

v_prop =  np.exp(-1j*V_x*dt/hbar_prim)
p_prop = np.exp(-1j*p**2*dt/(hbar_prim*2*m_h))

# set initial wavefunction
phi_x = wave_package_position(x, x0, p0, d, hbar_prim)
phi_x_initial = phi_x
n_x_initial = np.abs(phi_x_initial)**2


create_figures = False
fig, ax = plt.subplots() 
# propagate wave!
for t in range(Nt):
    phi_x = propagate(phi_x, p_prop, v_prop)
    n_x = np.abs(phi_x)**2
    
    if (create_figures and t%100 == 0):
        print(f'Saving figure for timestep {t} / {Nt}')
        ax.plot(x, n_x, color='blue')
        ax.plot(x, V_x, color='black')
        ax.set_xlim([-15, 15])
        ax.set_ylim([-0.01, 1.4])
        plt.savefig('../T3/figure/initial_gauss_momentum{:06d}.png'.format(t), bbox_inches='tight')
        ax.clear()

fig, ax = plt.subplots()
ax.plot(x, n_x_initial)
ax.plot(x, n_x)
    