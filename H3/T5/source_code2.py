import numpy as np
import matplotlib.pyplot as plt

## Create initial wave packet
def wave_package_position(x, x0, p0, d, hbar_prim):    
    prefactor = np.sqrt(np.sqrt(np.pi)*d)
    envelope  = np.exp(-(x - x0)**2/(2*d**2))
    plane_wave = np.exp(1j*p0*(x - x0)/hbar_prim)
    return (1/prefactor)*envelope*plane_wave

## Diabatic potential element V_11
def V11(x, a, b):
    output = np.zeros(len(x))
    output[x > 0] = a*(2-np.exp(-x[x > 0]/b))
    output[x <= 0] = a*np.exp(x[x <= 0]/b)
    return output

## Diabatic potential element V_22
def V22(x, a, b):
    return 2*a-V11(x, a, b)

## Diabatic potential element V_12
def V12(x, c, d, lim):
    v = c*np.exp(-(x/d)**2)
    v[v < lim] = lim
    return v

## Adiabatic potential element V_11. Ground state
def V11_adiabatic(x, a, b, c, d, lim):
    V_11 = V11(x, a, b)
    V_22 = V22(x, a, b)
    V_12 = V12(x, c, d, lim)
    return (V_11 + V_22 - np.sqrt((V_11-V_22)**2 + 4*V_12**2))/2 

## Adiabatic potential element V_22. Excited state
def V22_adiabatic(x, a, b, c, di, lim):
    V_11 = V11(x, a, b)
    V_22 = V22(x, a, b)
    V_12 = V12(x, c, d, lim)
    return (V_11 + V_22 + np.sqrt((V_11-V_22)**2 + 4*V_12**2))/2

## Function to propagate phi_x
def propagate(phi_x, p_prop, v_prop_1, v_prop_2, A, A_inv, Nx):
    ## Initialize matrices
    phi_trans_v = np.zeros_like(phi_x)
    phi_trans_t = np.zeros_like(phi_x)
    phi_ad = np.zeros_like(phi_x)
    phi_dia = np.zeros_like(phi_x)
    
    ## Transformation to adiabatic representation to apply potential operator
    phi_ad = np.einsum('ijk, ik -> jk', A, phi_x)
       
    ## Applying potential operator
    phi_trans_v[0,:] = v_prop_1*phi_ad[0,:]
    phi_trans_v[1,:] = v_prop_2*phi_ad[1,:]
    
    ## Transformation to diabatic representation to apply kinetic operator
    phi_dia = np.einsum('ijk, ik -> jk', A_inv, phi_trans_v)
    
    ## Applying kinetic operator
    phi_trans_t[0,:] = np.fft.ifft(p_prop*np.fft.fft(phi_dia[0,:]))
    phi_trans_t[1,:] = np.fft.ifft(p_prop*np.fft.fft(phi_dia[1,:]))
    return phi_trans_t

## Calculate the reflection
def reflection(phi_x, Nx, dx):
    n_x = np.abs(phi_x)**2
    return np.sum(n_x[0:int(Nx/2)-1])*dx

## Calculate the transmission
def transmission(phi_x, Nx, dx):
    n_x = np.abs(phi_x)**2
    return np.sum(n_x[int(Nx/2):-1])*dx

## Constants
hbar_si   = 1.054571817e-34           # hbar in SI
hbar_prim = hbar_si*1e34/1.602      # hbar in ASU_prim. eV*fs
m_prim_u  = 0.009647                 # mass/u in ASU_prim
m_h       = 1/m_prim_u                    # Mass of our hydrogen atom

## Grid
dx      = 0.01
x_start = -500.0
x_stop  = 500.0
x       = np.arange(x_start, x_stop, dx)
Nx      = len(x)

## Generating momentums for momentum operator
dp        = 2*np.pi*hbar_prim/(Nx*dx)              # Delta momentum
p_nyquist = 2*np.pi*hbar_prim/(2*dx)               # Nyquivst momentum
p         = np.arange(-p_nyquist, p_nyquist, dp)   # Generating momentum "freq"
p         = np.fft.fftshift(p)

## Time
T  = 2000
dt = 1
Nt = int(T/dt)

## Define kinetic operator
P_prop = np.exp(-1j*p**2*dt/(hbar_prim*2*m_h))

## Define potential energies
a = 0.3
b = 0.4
c = 0.1
d = 0.7 
lim = 1e-100    # Lower limit set for V_12 to not generate NaN

## Adiabatic potentials
V_a = V11_adiabatic(x, a, b, c, d, lim)     # V_a ground potental.
V_b = V22_adiabatic(x, a, b, c, d, lim)     # V_b excited potential.

## Define potential operators
V_prop_1 = np.exp(-1j*V_a*dt/hbar_prim)
V_prop_2 = np.exp(-1j*V_b*dt/hbar_prim)

## Diabatic potentials
V_11 = V11(x, a, b)
V_22 = V22(x, a, b)
V_12 = V12(x, c, d, lim)
V_dia = np.array([[V_11, V_12], [V_12, V_22]])

# Create transformation matrix
alpha = V_22 - V_11
beta = np.sqrt((V_11-V_22)**2 + 4*V_12**2)

A_11 = -(alpha+beta)/(2*V_12)
A_12 = -(alpha-beta)/(2*V_12)
A_21 = np.ones(len(x))
A_22 = np.ones(len(x))
A = np.array([[A_11, A_12], [A_21, A_22]])
A_inv = -(V_12/beta)*np.array([[A_22, -A_12], [-A_21, A_11]])

## Check that the matrix A diagonalizes V_dia
V_ad_check = np.array([[V_a, np.zeros(len(x))], [np.zeros(len(x)), V_b]])
V_ad = np.zeros_like(V_dia)
for i in range(len(x)):
    V_ad[:,:,i] = (A_inv[:,:,i] @ V_dia[:,:,i] @ A[:,:,i])

print(np.allclose(V_ad, V_ad_check, atol=1e-15))

## Initial values for wavefunction
E_i = 0.25                                     # initial energy
width = hbar_prim*np.sqrt(1/(0.04*E_i*m_h))    # Width of our hydrogen atom
p0 = np.sqrt(2*m_h*E_i)     # Initial momentum of our hydrogen atom
x0 = -15                    # Initial position of our hydrogen atom

## Set initial wavefunction
phi_x_1 = wave_package_position(x, x0, p0, width, hbar_prim)
phi_x_2 = np.zeros(len(x))
phi_x = np.array([phi_x_1, phi_x_2])

## Propagate system!
for t in range(Nt):
    if (t%100 == 0):
        print(f'Timestep {t} \ {Nt}') 
    phi_x = propagate(phi_x, P_prop, V_prop_1, V_prop_2, A, A_inv, Nx)
       
## Calculate final transmission and refection     
T_1 = transmission(phi_x[0,:], Nx, dx)
T_2 = transmission(phi_x[1,:], Nx, dx)
R_1 = reflection(phi_x[0,:], Nx, dx)
R_2 = reflection(phi_x[1,:], Nx, dx)
