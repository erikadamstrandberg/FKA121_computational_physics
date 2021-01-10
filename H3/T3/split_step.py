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
x_start = -150.0
x_stop  = 150.0
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

initial_energy = 0.08
p0 = np.sqrt(2*initial_energy*m_h)   # Initial momentum of our hydrogen atom
x0 = -15                 # Initial position of our hydrogen atom

## Propagation
T = 1e6
dt = 0.5
Nt = int(T/dt)

V0 = 0.1
alpha = 2.0

# adiabatic potentials
lim = 1e-60
a = 0.3
b = 0.4
c = 0.05

V_x = V0/(np.cosh(x/alpha)**2)

v_prop =  np.exp(-1j*V_x*dt/hbar_prim)
p_prop = np.exp(-1j*p**2*dt/(hbar_prim*2*m_h))

# set initial wavefunction
phi_x = wave_package_position(x, x0, p0, d, hbar_prim)
phi_x_initial = phi_x
n_x_initial = np.abs(phi_x_initial)**2

fig, ax = plt.subplots() 
ax.plot(x, V_x        , color='orange', linewidth=2, label=r'$V_x(x)$')
ax.plot(x, n_x_initial, color='blue', linewidth=2, label=r'$|\psi(x)|^2$')

ax.set_title(r'Simulation set up, $V_0=0.1$ eV, $\alpha=2.0$ Å', fontsize='16')
ax.set_xlabel(r'$x$ [$Å$]', fontsize='16')
ax.set_ylabel(r'$|\psi(x)|^2$ [$Å^{-1}$], $V(x)$ [$eV$]', fontsize='16')

ax.set_xlim([-20,20])

ax.grid()
ax.legend(fontsize='16', loc='upper right')

#plt.savefig('../T3/Sim_set_up_alpha_2.pdf', format='pdf', bbox_inches='tight')
plt.show()
 
V_less_then = 1e-8
V_stop_index = np.argmax(V_x > V_less_then)

x_stop_index_R = np.argmax(x > -100)
x_stop = x[V_stop_index]
x_stop_2 = x[-1-V_stop_index]
y_stop = np.array([0,0.3])

#%%
stop_prop = 1e-6
start_checking = False
# propagate wave!
for t in range(Nt):
    if (t%100 == 0):
        print(f'Saving figure for timestep: {t} / {Nt}')
        print(np.abs(phi_x[x_stop_index_R])**2)
    if (np.abs(phi_x[x_stop_index_R])**2 > stop_prop):
        break
    phi_x = propagate(phi_x, p_prop, v_prop)
    
n_x = np.abs(phi_x)**2

fig, ax = plt.subplots() 
ax.plot(x, V_x        , color='orange', linewidth=2, label=r'$V_x(x)$')
ax.plot(x, n_x, color='blue', linewidth=2, label=r'$|\psi(x)|^2$')

print(n_x[x_stop_index_R])
print(n_x[V_stop_index])
print(n_x[-1-V_stop_index])

#np.savetxt('data/V_x_e0_{:.2f}_V0_{:.2f}.csv'.format(initial_energy,V0), V_x, delimiter=",")
#np.savetxt('data/n_x_e0_{:.2f}_V0_{:.2f}.csv'.format(initial_energy,V0), n_x, delimiter=",")
#np.savetxt('data/x_e0_{:.2f}_V0_{:.2f}.csv'.format(initial_energy,V0), x, delimiter=",")

    