#%%
import numpy as np
import matplotlib.pyplot as plt

def wave_package_position(x, x0, p0, d, hbar_prim):    
    prefactor = np.sqrt(np.sqrt(np.pi)*d)
    envelope  = np.exp(-(x - x0)**2/(2*d**2))
    plane_wave = np.exp(1j*p0*(x - x0)/hbar_prim)
    return (1/prefactor)*envelope*plane_wave

def propagate(phi_x, p_prop, v_prop):
    return np.fft.ifft(p_prop*np.fft.fft(v_prop*phi_x))

def reflection(n_x, Nx, dx):
    return np.sum(n_x[0:int(Nx/2)-1])*dx

def transmission(n_x, Nx, dx):
    return np.sum(n_x[int(Nx/2):-1])*dx


#%%
    
x = np.loadtxt('../T3/data/x_e0_0.08_V0_0.10.csv', delimiter=',')
n_x = np.loadtxt('../T3/data/n_x_e0_0.08_V0_0.10.csv', delimiter=',')
V_x = np.loadtxt('../T3/data/V_x_e0_0.08_V0_0.10.csv', delimiter=',')

Nx = len(x)
dx = x[1]-x[0]

fig, ax = plt.subplots() 
ax.plot(x, V_x, color='orange', linewidth=2, label=r'$V(x)$')
ax.plot(x, n_x, color='blue', linewidth=2, label=r'$|\psi(x)|^2$')

T = transmission(n_x, Nx, dx)
R = reflection(n_x, Nx, dx)

print(f'Transsmission: {T}')
print(f'Reflection:    {R}')

ax.set_title(r'After simulation, $\alpha=0.5$ $Å$, $p_0^2/2m = 0.08$ $eV$', fontsize='16')
ax.set_xlabel(r'$x$ [$Å$]', fontsize='16')
ax.set_ylabel(r'$|\psi(x)|^2$ [$Å^{-1}$], $V(x)$ [eV]', fontsize='16')

ax.grid()
ax.set_xlim([-200,200])
ax.legend(fontsize='16', loc='upper right')

#%%
    
x = np.loadtxt('../T3/data/x_e0_0.10_V0_0.10.csv', delimiter=',')
n_x = np.loadtxt('../T3/data/n_x_e0_0.10_V0_0.10.csv', delimiter=',')
V_x = np.loadtxt('../T3/data/V_x_e0_0.10_V0_0.10.csv', delimiter=',')

Nx = len(x)
dx = x[1]-x[0]

fig, ax = plt.subplots() 
ax.plot(x, V_x, color='orange', linewidth=2, label=r'$V(x)$')
ax.plot(x, n_x, color='blue', linewidth=2, label=r'$|\psi(x)|^2$')

ax.set_xlim([-200,200])

T = transmission(n_x, Nx, dx)
R = reflection(n_x, Nx, dx)

print(f'Transsmission: {T}')
print(f'Reflection:    {R}')

ax.set_title(r'After simulation, $\alpha=0.5$ $Å$, $p_0^2/2m = 0.10$ $eV$', fontsize='16')
ax.set_xlabel(r'$x$ [$Å$]', fontsize='16')
ax.set_ylabel(r'$|\psi(x)|^2$ [$Å^{-1}$], $V(x)$ [eV]', fontsize='16')

ax.grid()
ax.set_xlim([-200,200])
ax.legend(fontsize='16', loc='upper right')

#%%
    
x = np.loadtxt('../T3/data/x_e0_0.12_V0_0.10.csv', delimiter=',')
n_x = np.loadtxt('../T3/data/n_x_e0_0.12_V0_0.10.csv', delimiter=',')
V_x = np.loadtxt('../T3/data/V_x_e0_0.12_V0_0.10.csv', delimiter=',')

Nx = len(x)
dx = x[1]-x[0]

fig, ax = plt.subplots() 
ax.plot(x, V_x, color='orange', linewidth=2, label=r'$V(x)$')
ax.plot(x, n_x, color='blue', linewidth=2, label=r'$|\psi(x)|^2$')

T = transmission(n_x, Nx, dx)
R = reflection(n_x, Nx, dx)

print(f'Transsmission: {T}')
print(f'Reflection:    {R}')

ax.set_title(r'After simulation, $\alpha=0.5$ $Å$, $p_0^2/2m = 0.12$ $eV$', fontsize='16')
ax.set_xlabel(r'$x$ [$Å$]', fontsize='16')
ax.set_ylabel(r'$|\psi(x)|^2$ [$Å^{-1}$], $V(x)$ [eV]', fontsize='16')

ax.grid()
ax.set_xlim([-200,200])
ax.legend(fontsize='16', loc='upper right')

#%%
    
x = np.loadtxt('../T3/data/x_e0_0.08_V0_0.10_wide.csv', delimiter=',')
n_x = np.loadtxt('../T3/data/n_x_e0_0.08_V0_0.10_wide.csv', delimiter=',')
V_x = np.loadtxt('../T3/data/V_x_e0_0.08_V0_0.10_wide.csv', delimiter=',')

Nx = len(x)
dx = x[1]-x[0]

fig, ax = plt.subplots() 
ax.plot(x, V_x, color='orange', linewidth=2, label=r'$V(x)$')
ax.plot(x, n_x, color='blue', linewidth=2, label=r'$|\psi(x)|^2$')

T = transmission(n_x, Nx, dx)
R = reflection(n_x, Nx, dx)

print(f'Transsmission: {T}')
print(f'Reflection:    {R}')

ax.set_title(r'After simulation, $\alpha=2.0$ $Å$, $p_0^2/2m = 0.08$ $eV$', fontsize='16')
ax.set_xlabel(r'$x$ [$Å$]', fontsize='16')
ax.set_ylabel(r'$|\psi(x)|^2$ [$Å^{-1}$], $V(x)$ [eV]', fontsize='16')

ax.grid()
ax.set_xlim([-120,120])
ax.legend(fontsize='16', loc='upper right')



#%%
    
x = np.loadtxt('../T3/data/x_e0_0.10_V0_0.10_wide.csv', delimiter=',')
n_x = np.loadtxt('../T3/data/n_x_e0_0.10_V0_0.10_wide.csv', delimiter=',')
V_x = np.loadtxt('../T3/data/V_x_e0_0.10_V0_0.10_wide.csv', delimiter=',')

Nx = len(x)
dx = x[1]-x[0]

fig, ax = plt.subplots() 
ax.plot(x, V_x, color='orange', linewidth=2, label=r'$V(x)$')
ax.plot(x, n_x, color='blue', linewidth=2, label=r'$|\psi(x)|^2$')

T = transmission(n_x, Nx, dx)
R = reflection(n_x, Nx, dx)

print(f'Transsmission: {T}')
print(f'Reflection:    {R}')

ax.set_title(r'After simulation, $\alpha=2.0$ $Å$, $p_0^2/2m = 0.10$ $eV$', fontsize='16')
ax.set_xlabel(r'$x$ [$Å$]', fontsize='16')
ax.set_ylabel(r'$|\psi(x)|^2$ [$Å^{-1}$], $V(x)$ [eV]', fontsize='16')

ax.grid()
ax.set_xlim([-120,120])
ax.legend(fontsize='16', loc='upper right')


#%%
    
x = np.loadtxt('../T3/data/x_e0_0.12_V0_0.10_wide.csv', delimiter=',')
n_x = np.loadtxt('../T3/data/n_x_e0_0.12_V0_0.10_wide.csv', delimiter=',')
V_x = np.loadtxt('../T3/data/V_x_e0_0.12_V0_0.10_wide.csv', delimiter=',')

Nx = len(x)
dx = x[1]-x[0]

fig, ax = plt.subplots() 
ax.plot(x, V_x, color='orange', linewidth=2, label=r'$V(x)$')
ax.plot(x, n_x, color='blue', linewidth=2, label=r'$|\psi(x)|^2$')

T = transmission(n_x, Nx, dx)
R = reflection(n_x, Nx, dx)

print(f'Transsmission: {T}')
print(f'Reflection:    {R}')

ax.set_title(r'After simulation, $\alpha=2.0$ $Å$, $p_0^2/2m = 0.12$ $eV$', fontsize='16')
ax.set_xlabel(r'$x$ [$Å$]', fontsize='16')
ax.set_ylabel(r'$|\psi(x)|^2$ [$Å^{-1}$], $V(x)$ [eV]', fontsize='16')

ax.grid()
ax.set_xlim([-120,120])
ax.legend(fontsize='16', loc='upper right')


