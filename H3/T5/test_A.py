import numpy as np
import matplotlib.pyplot as plt
import time

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

def V22_adiabatic(x, a, b, c, d, lim):
    V_11 = V11(x, a, b)
    V_22 = V22(x, a, b)
    V_12 = V12(x, c, d, lim)
    return (V_11 + V_22 - np.sqrt((V_11-V_22)**2 + 4*V_12**2))/2 


#%% define constants
a = 0.3
b = 0.4
c = 0.05
d = 0.7

dx      = 0.001
x_start = -5.0
x_stop  = 5.0
x = np.arange(x_start, x_stop, dx)

lim = 1e-60
# diabatic potentials
V_11 = V11(x, a, b)
V_22 = V22(x, a, b)
V_12 = V12(x, c, d, lim)

# adiabatic potentials for coupling c1=0.05 and c2=0.1
V_a = V11_adiabatic(x, a, b, c, d, lim)
V_b = V22_adiabatic(x, a, b, c, d, lim)

# transformation matrix
alpha = V_22 - V_11
beta = np.sqrt((V_11-V_22)**2 + 4*V_12**2)

A_11 = -(alpha-beta)/(2*V_12)
A_12 = -(alpha+beta)/(2*V_12)
A_21 = np.ones(len(x))
A_22 = np.ones(len(x))
A = np.array([[A_11, A_12], [A_21, A_22]])
A_inv = (V_12/beta)*np.array([[A_22, -A_12], [-A_21, A_11]])

# transform from diabatic to adiabatic
V_dia = np.array([[V_11, V_12], [V_12, V_22]])

start = time.time()
V_ad = np.einsum('ijn,jkn->ikn', A_inv, np.einsum('ijn,jkn->ikn', V_dia, A))
stop = time.time()
print(stop-start)

V_ad = np.zeros_like(V_dia)
start_2 = time.time()
for i in range(len(x)):
    V_ad[:,:,i] = A_inv[:,:,i] @ V_dia[:,:,i] @ A[:,:,i]
stop_2 = time.time()
print(stop_2-start_2)

plt.plot(x, V_a)
plt.plot(x, V_ad[0,0,:])
plt.show()

