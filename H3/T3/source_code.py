import numpy as np
import matplotlib.pyplot as plt

def wave_package_position(x, x0, p0, d, hbar_prim):    
    prefactor = np.sqrt(np.sqrt(np.pi)*d)
    envelope  = np.exp(-(x - x0)**2/(2*d**2))
    plane_wave = np.exp(1j*p0*(x - x0)/hbar_prim)
    return (1/prefactor)*envelope*plane_wave

def propagate(phi_x, p_prop, v_prop):
    return np.fft.ifft(p_prop*np.fft.fft(v_prop*phi_x))

## Constants
hbar_si   = 1.054571817e-34           # hbar in SI
hbar_prim = hbar_si*1e34/1.602        # hbar in ASU_prim. eV*fs
m_prim_u  = 0.009647                  # mass/u in ASU_prim

## Grid
dx      = 0.001         # dt
x_start = -150.0
x_stop  = 150.0
x       = np.arange(x_start, x_stop, dx)
Nx      = len(x)

## Generating momentums for momentum operator
dp        = 2*np.pi*hbar_prim/(Nx*dx)              # Delta momentum
p_nyquist = 2*np.pi*hbar_prim/(2*dx)               # Nyquivst momentum
p         = np.arange(-p_nyquist, p_nyquist, dp)   # Generating momentum "freq"
p         = np.fft.fftshift(p)

## Initial values for wavefunction
d   = 0.5               # Width of our hydrogen atom
m_h = 1/m_prim_u        # Mass of our hydrogen atom
initial_energy = 0.08
p0 = np.sqrt(2*initial_energy*m_h)   # Initial momentum of our hydrogen atom
x0 = -15                             # Initial position of our hydrogen atom

## Time
T = 1e6             # Total time we allow simulation to run for
dt = 0.5            # dt
Nt = int(T/dt)

## Eckart barrier 
V0 = 0.1
alpha = 0.5
V_x = V0/(np.cosh(x/alpha)**2)

## Operators
v_prop =  np.exp(-1j*V_x*dt/hbar_prim)
p_prop = np.exp(-1j*p**2*dt/(hbar_prim*2*m_h))

# Generate initial wavefunction
phi_x = wave_package_position(x, x0, p0, d, hbar_prim)
phi_x_initial = phi_x
n_x_initial = np.abs(phi_x_initial)**2

## Set stopping condition for simulation
## When the reflected proababilty reached over 1e-6 at x=-200 the simulation is stopped
x_stop_index_R = np.argmax(x > -150)
stop_prop = 1e-6


## Propagate system
for t in range(Nt):
    if (t%100 == 0):
        print(f'Saving figure for timestep: {t} / {Nt}')
    ## Stops simulation
    if (np.abs(phi_x[x_stop_index_R])**2 > stop_prop):
        break
    phi_x = propagate(phi_x, p_prop, v_prop)
    
n_x = np.abs(phi_x)**2
