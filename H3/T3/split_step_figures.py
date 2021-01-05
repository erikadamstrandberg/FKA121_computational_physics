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
    
x = np.loadtxt('data/x_e0_0.10_V0_0.10.csv', delimiter=',')
n_x = np.loadtxt('data/n_x_e0_0.10_V0_0.10.csv', delimiter=',')
V_x = np.loadtxt('data/V_x_e0_0.10_V0_0.10.csv', delimiter=',')

Nx = len(x)
dx = x[1]-x[0]

fig, ax = plt.subplots() 
ax.plot(x, n_x)
ax.plot(x, V_x)

ax.set_xlim([-100,100])

T = transmission(n_x, Nx, dx)
R = reflection(n_x, Nx, dx)

print(f'Transsmission: {T}')
print(f'Reflection:    {R}')


#%%
    
x = np.loadtxt('data/x_e0_0.08_V0_0.10.csv', delimiter=',')
n_x = np.loadtxt('data/n_x_e0_0.08_V0_0.10.csv', delimiter=',')
V_x = np.loadtxt('data/V_x_e0_0.08_V0_0.10.csv', delimiter=',')

Nx = len(x)
dx = x[1]-x[0]

fig, ax = plt.subplots() 
ax.plot(x, n_x)
ax.plot(x, V_x)

ax.set_xlim([-100,100])

T = transmission(n_x, Nx, dx)
R = reflection(n_x, Nx, dx)

print(f'Transsmission: {T}')
print(f'Reflection:    {R}')

#%%
    
x = np.loadtxt('data/x_e0_0.12_V0_0.10.csv', delimiter=',')
n_x = np.loadtxt('data/n_x_e0_0.12_V0_0.10.csv', delimiter=',')
V_x = np.loadtxt('data/V_x_e0_0.12_V0_0.10.csv', delimiter=',')

Nx = len(x)
dx = x[1]-x[0]

fig, ax = plt.subplots() 
ax.plot(x, n_x)
ax.plot(x, V_x)

ax.set_xlim([-100,100])

T = transmission(n_x, Nx, dx)
R = reflection(n_x, Nx, dx)

print(f'Transsmission: {T}')
print(f'Reflection:    {R}')


#%%
    
x = np.loadtxt('data/x_e0_0.08_V0_0.10.csv', delimiter=',')
n_x = np.loadtxt('data/n_x_e0_0.12_V0_0.10.csv', delimiter=',')
V_x = np.loadtxt('data/V_x_e0_0.12_V0_0.10.csv', delimiter=',')

Nx = len(x)
dx = x[1]-x[0]

fig, ax = plt.subplots() 
ax.plot(x, n_x)
ax.plot(x, V_x)

ax.set_xlim([-100,100])

T = transmission(n_x, Nx, dx)
R = reflection(n_x, Nx, dx)

print(f'Transsmission: {T}')
print(f'Reflection:    {R}')


