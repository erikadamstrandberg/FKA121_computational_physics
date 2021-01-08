import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv

def wave_package_position(x, x0, p0, d, hbar_prim):    
    prefactor = np.sqrt(np.sqrt(np.pi)*d)
    envelope  = np.exp(-(x - x0)**2/(2*d**2))
    plane_wave = np.exp(1j*p0*(x - x0)/hbar_prim)
    return (1/prefactor)*envelope*plane_wave

def V11(x, a, b):
    output = np.zeros(len(x))
    output[x > 0] = a*(2-np.exp(-x[x > 0]/b))
    output[x <= 0] = a*np.exp(x[x <= 0]/b)
    return output

def V22(x, a, b):
    return 2*a-V11(x, a, b)

def V12(x, c, d, lim):
    v = c*np.exp(-(x/d)**2)
    v[v < lim] = lim
    return v

def V11_adiabatic(x, a, b, c, d, lim):
    V_11 = V11(x, a, b)
    V_22 = V22(x, a, b)
    V_12 = V12(x, c, d, lim)
    return (V_11 + V_22 + np.sqrt((V_11-V_22)**2 + 4*V_12**2))/2 

def V22_adiabatic(x, a, b, c, di, lim):
    V_11 = V11(x, a, b)
    V_22 = V22(x, a, b)
    V_12 = V12(x, c, d, lim)
    return (V_11 + V_22 - np.sqrt((V_11-V_22)**2 + 4*V_12**2))/2

####################################################################

def linalg_inv(A):
    A_inv = np.zeros_like(A)
    for i in range(len(A[0,0,:])):
        A_inv[:,:,i] = inv(A[:,:,i])
    return A_inv

def to_diabatic(phi_x, A):
    return np.einsum('ijk,ik->jk', A, phi_x)

def to_adiabatic(phi_x, A_inv):
    return np.einsum('ijk,ik->jk', A_inv, phi_x)

def propagate(phi_x, p_prop, v_prop_1, v_prop_2, A, A_inv, Nx):
    phi_trans_v = np.zeros_like(phi_x)
    phi_trans_t = np.zeros_like(phi_x)
    
    phi_trans_v[0,:] = v_prop_1*phi_x[0,:]
    phi_trans_v[1,:] = v_prop_2*phi_x[1,:]
    
    #phi_dia = np.zeros_like(phi_x)
    #phi_dia = np.einsum('ijk,ik->jk', A, phi_trans_v)
    phi_dia = np.zeros_like(phi_x)
    for i in range(Nx):
        phi_dia[:,i] = A[:,:,i] @ phi_trans_v[:,i]
        
    phi_trans_t[0,:] = np.fft.ifft(p_prop*np.fft.fft(phi_dia[0,:]))
    phi_trans_t[1,:] = np.fft.ifft(p_prop*np.fft.fft(phi_dia[1,:]))
    
    phi_ad = np.zeros_like(phi_x)
    for i in range(Nx):
        phi_ad[:,i] = A_inv[:,:,i] @ phi_trans_t[:,i]
    #phi_ad =  np.einsum('ijk,ik->jk', A_inv, phi_trans_t)
    return phi_ad # np.einsum('ijk,ik->jk', A_inv, phi_trans_t)
    

#%% define constants

hbar_si = 1.054571817e-34           # hbar in SI
hbar_prim = hbar_si*1e34/1.602      # hbar in ASU_prim. eV*fs
m_prim_u = 0.009647                 # mass/u in ASU_prim


# position space
dx      = 0.01
x_start = -10.0
x_stop  = 10.0
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

p0 = np.sqrt(0.2*m_h)   # Initial momentum of our hydrogen atom
x0 = -7                 # Initial position of our hydrogen atom

#%% Propagation
T = 120
dt = 0.1
Nt = int(T/dt)


# adiabatic potentials
lim = 1e-60
a = 0.3
b = 0.4
c = 0.05


V_a = V11_adiabatic(x, a, b, c, d, lim) 
V_b = V22_adiabatic(x, a, b, c, d, lim) 

V_prop_1 =  np.exp(-1j*V_a*dt/hbar_prim)
V_prop_2 =  np.exp(-1j*V_b*dt/hbar_prim)

P_prop = np.exp(-1j*p**2*dt/(hbar_prim*2*m_h))

# create transformation matrix
# diabatic potentials
V_11 = V11(x, a, b)
V_22 = V22(x, a, b)
V_12 = V12(x, c, d, lim)

alpha = V_22 - V_11
beta = np.sqrt((V_11-V_22)**2 + 4*V_12**2)

A_11 = -(alpha-beta)/(2*V_12)
A_12 = -(alpha+beta)/(2*V_12)
A_21 = np.ones(len(x))
A_22 = np.ones(len(x))
A = np.array([[A_11, A_12], [A_21, A_22]])

A_inv = (V_12/beta)*np.array([[A_22, -A_12], [-A_21, A_11]])
#A_inv = linalg_inv(A)
# set initial wavefunction
phi_x_1 = wave_package_position(x, x0, p0, d, hbar_prim)
phi_x_2 = np.zeros(len(x))
phi_x = np.array([phi_x_1, phi_x_2])

phi_x_initial_1 = phi_x[0,:]
#phi_x_initial_2 = phi_x[1,:]


fig, ax = plt.subplots()

# propagate wave!
#phi_x_prop = propagate(phi_x, P_prop, V_prop_1, V_prop_2, A, A_inv)
for t in range(Nt):
    phi_x = propagate(phi_x, P_prop, V_prop_1, V_prop_2, A, A_inv, Nx)


ax.plot(x, phi_x[0,:])
ax.plot(x, phi_x[1,:])
#ax.plot(x, phi_x_initial_1)

plt.show()
