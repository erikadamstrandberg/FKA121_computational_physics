import numpy as np
import matplotlib.pyplot as plt

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
    return (V_11 + V_22 - np.sqrt((V_11-V_22)**2 + 4*V_12**2))/2 

def V22_adiabatic(x, a, b, c, di, lim):
    V_11 = V11(x, a, b)
    V_22 = V22(x, a, b)
    V_12 = V12(x, c, d, lim)
    return (V_11 + V_22 + np.sqrt((V_11-V_22)**2 + 4*V_12**2))/2

def propagate(phi_x, p_prop, v_prop_1, v_prop_2, A, A_inv, Nx):
    phi_trans_v = np.zeros_like(phi_x)
    phi_trans_t = np.zeros_like(phi_x)
    phi_dia = np.zeros_like(phi_x)
    
    # Propagating with potential operatior
    phi_trans_v[0,:] = v_prop_1*phi_x[0,:]
    phi_trans_v[1,:] = v_prop_2*phi_x[1,:]
    
    # Transformation to diabatic states
    phi_dia = np.einsum('ijk, ik -> jk', A_inv, phi_trans_v)
    
    # Propagating with kinetic operatior
    phi_trans_t[0,:] = np.fft.ifft(p_prop*np.fft.fft(phi_dia[0,:]))
    phi_trans_t[1,:] = np.fft.ifft(p_prop*np.fft.fft(phi_dia[1,:]))
    
    # Transformation back to adiabatic states
    return np.einsum('ijk, ik -> jk', A, phi_trans_t)

def reflection(phi_x, Nx, dx):
    n_x = np.abs(phi_x)**2
    return np.sum(n_x[0:int(Nx/2)-1])*dx

def transmission(phi_x, Nx, dx):
    n_x = np.abs(phi_x)**2
    return np.sum(n_x[int(Nx/2):-1])*dx

    

#%% define constants

hbar_si = 1.054571817e-34           # hbar in SI
hbar_prim = hbar_si*1e34/1.602      # hbar in ASU_prim. eV*fs
m_prim_u = 0.009647                 # mass/u in ASU_prim
m_h = 1/m_prim_u                    # Mass of our hydrogen atom

# position space
dx      = 0.01
x_start = -80.0
x_stop  = 80.0
x = np.arange(x_start, x_stop, dx)
Nx = len(x)

# momentum space
dp        = 2*np.pi*hbar_prim/(Nx*dx)              # Delta momentum
p_nyquist = 2*np.pi*hbar_prim/(2*dx)               # Nyquivst momentum
p         = np.arange(-p_nyquist, p_nyquist, dp)   # Generating momentum "freq"
p         = np.fft.fftshift(p)

# time propagation
T = 300
dt = 0.1
Nt = int(T/dt)

# Define potential energies
lim = 1e-100
a = 0.3
b = 0.4
c = 0.05
d = 0.7 

# adiabatic potentials
V_a = V11_adiabatic(x, a, b, c, d, lim)     # V_a ground potental.
V_b = V22_adiabatic(x, a, b, c, d, lim)     # V_b excited potential.

# diabatic potentials
V_11 = V11(x, a, b)
V_22 = V22(x, a, b)
V_12 = V12(x, c, d, lim)

# create transformation matrix
alpha = V_22 - V_11
beta = np.sqrt((V_11-V_22)**2 + 4*V_12**2)

A_11 = -(alpha+beta)/(2*V_12)
A_12 = -(alpha-beta)/(2*V_12)
A_21 = np.ones(len(x))
A_22 = np.ones(len(x))
A = V_12*np.array([[A_11, A_12], [A_21, A_22]])
A_inv = -(1/beta)*np.array([[A_22, -A_12], [-A_21, A_11]])


V_ad_check = np.array([[V_a, np.zeros(len(x))], [np.zeros(len(x)), V_b]])
V_dia = np.array([[V_11, V_12], [V_12, V_22]])

V_ad = np.zeros_like(V_dia)
for i in range(len(x)):
    V_ad[:,:,i] = (A_inv[:,:,i] @ V_dia[:,:,i] @ A[:,:,i])

print(np.allclose(V_ad, V_ad_check, atol=1e-15))

# define propagators
V_prop_1 =  np.exp(-1j*V_a*dt/hbar_prim)
V_prop_2 =  np.exp(-1j*V_b*dt/hbar_prim)

P_prop = np.exp(-1j*p**2*dt/(hbar_prim*2*m_h))

# initial values for wavefunction
E_i = 0.6*a                                     # initial energy
width = hbar_prim/(4*np.sqrt(0.1*E_i*m_h))    # Width of our hydrogen atom

p0 = np.sqrt(2*m_h*E_i)     # Initial momentum of our hydrogen atom
x0 = -15                    # Initial position of our hydrogen atom

# set initial wavefunction
phi_x_1 = wave_package_position(x, x0, p0, width, hbar_prim)
phi_x_2 = np.zeros(len(x))
phi_x = np.array([phi_x_1, phi_x_2])

T_1 = transmission(phi_x[0,:], Nx, dx)
T_2 = transmission(phi_x[1,:], Nx, dx)
R_1 = reflection(phi_x[0,:], Nx, dx)
R_2 = reflection(phi_x[1,:], Nx, dx)
prob_0 = T_1 + T_2 + R_1 + R_2
print(prob_0)

x_stop_R1 = np.argmax(x > -200)
V_check_index = np.argmax(V_b > 1e-8)
wave_packet_check = 1e-10

# propagate wave!
for t in range(Nt):
    if (t%100 == 0):
        print(f'Timestep {t} \ {Nt}') 
        if (np.abs(phi_x[0,x_stop_R1])**2 > wave_packet_check):
            break
    phi_x = propagate(phi_x, P_prop, V_prop_1, V_prop_2, A, A_inv, Nx)
    
    
T_1 = transmission(phi_x[0,:], Nx, dx)
T_2 = transmission(phi_x[1,:], Nx, dx)
R_1 = reflection(phi_x[0,:], Nx, dx)
R_2 = reflection(phi_x[1,:], Nx, dx)

print('T_1 = ' + str(T_1))
print('R_1 = ' + str(R_1))
print('T_2 = ' + str(T_2))
print('R_2 = ' + str(R_2))
print('total probability = ' + str(T_1 + T_2 + R_1 + R_2))

fig, ax = plt.subplots()
ax.plot(x, np.abs(phi_x[0,:])**2, label=r'$\phi_1(x)$')
ax.plot(x, np.abs(phi_x[1,:])**2, label=r'$\phi_2(x)$')
ax.legend()

plt.show()
