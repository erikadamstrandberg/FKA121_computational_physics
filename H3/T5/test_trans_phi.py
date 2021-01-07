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

#%% define constants

hbar_si = 1.054571817e-34           # hbar in SI
hbar_prim = hbar_si*1e34/1.602      # hbar in ASU_prim. eV*fs
m_prim_u = 0.009647                 # mass/u in ASU_prim


# position space
dx      = 0.001
x_start = -20.0
x_stop  = 20.0
x = np.arange(x_start, x_stop, dx)
Nx = len(x)

# initial values for wavefunction
d   = 0.5               # Width of our hydrogen atom
m_h = 1/m_prim_u        # Mass of our hydrogen atom

p0 = np.sqrt(0.2*m_h)   # Initial momentum of our hydrogen atom
x0 = -1                 # Initial position of our hydrogen atom

# adiabatic potentials
lim = 1e-10
a = 0.3
b = 0.4
c = 0.05

V_a = V11_adiabatic(x, a, b, c, d, lim) 
V_b = V22_adiabatic(x, a, b, c, d, lim) 

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

# set initial wavefunction
phi_x_1 = wave_package_position(x, x0, p0, d, hbar_prim)
phi_x_2 = np.zeros(len(x))
phi_x = np.array([[phi_x_1, phi_x_2]])

phi_trans = np.zeros_like(phi_x)
phi_back = np.zeros_like(phi_x)     # should have same values as phi_x

# transform from adiabatic to diabatic
for i in range(len(x)):
    phi_temp = phi_x[:,:,i] # (1,2)

    A_temp = A[:,:,0]   # (2,2)

    phi_trans_temp = phi_temp @ A_temp  # (1,2)
    phi_trans[:,:,i] = phi_trans_temp

# transform from diabatic to adiabatic
for i in range(len(x)):
    phi_temp = phi_trans[:,:,i] # (1,2)

    A_temp = A[:,:,0]   # (2,2)
    A_inv_temp = inv(A_temp) #A_inv[:,:,i]   # (2,2)

    phi_back_temp = phi_temp @ A_inv_temp # (1,2) 
    phi_back[:,:,i] = phi_back_temp

#print(phi_x)
#print(phi_trans)
#print(phi_back)


fig, ax = plt.subplots(1,2)
ax[0].plot(x, phi_x[0,0,:], label=r'$\phi_1^{ad}(x)$')
ax[0].plot(x, phi_x[0,1,:], label=r'$\phi_2^{ad}(x)$')
ax[0].plot(x, phi_back[0,0,:], label=r'$\phi_1^{back}$')
ax[0].plot(x, phi_back[0,1,:], label=r'$\phi_2^{back}$')
ax[0].legend()

ax[1].plot(x, phi_trans[0,0,:], label=r'$\phi_1^{di}(x)$')
ax[1].plot(x, phi_trans[0,1,:], label=r'$\phi_2^{di}(x)$')
ax[1].legend()
plt.show()


"""
##### test diagonalization of V #########################
# transform from diabatic to adiabatic
V_dia = np.array([[V_11, V_12], [V_12, V_22]])

# using einsum
V_ad_1 = np.einsum('ijn,jkn->ikn', A_inv, np.einsum('ijn,jkn->ikn', V_dia, A))
# using loop
V_ad_2 = np.zeros_like(V_dia)
for i in range(len(x)):
    V_ad_2[:,:,i] = A_inv[:,:,i] @ V_dia[:,:,i] @ A[:,:,i]

#plt.plot(x, V_b)
#plt.plot(x, V_ad_2[1,1,:])
#plt.show()
#########################################################
"""

