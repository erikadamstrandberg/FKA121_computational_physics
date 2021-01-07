import numpy as np
import matplotlib.pyplot as plt

def V11(x, a, b):
    output = np.zeros(len(x))
    output[x > 0] = a*(2-np.exp(-x[x > 0]/b))
    output[x <= 0] = a*np.exp(x[x <= 0]/b)
    return output

def V22(x, a, b):
    return 2*a-V11(x, a, b)

def V12(x, c, d):
    return c*np.exp(-(x/d)**2)

def V11_adiabatic(x, a, b, c, d):
    V_11 = V11(x, a, b)
    V_22 = V22(x, a, b)
    V_12 = V12(x, c, d)
    return (V_11 + V_22 + np.sqrt((V_11-V_22)**2 + 4*V_12**2))/2 

def V22_adiabatic(x, a, b, c, d):
    V_11 = V11(x, a, b)
    V_22 = V22(x, a, b)
    V_12 = V12(x, c, d)
    return (V_11 + V_22 - np.sqrt((V_11-V_22)**2 + 4*V_12**2))/2 


#%% define constants
a = 0.3
b = 0.4
c1 = 0.05
c2 = 0.1
d = 0.7

dx      = 0.001
x_start = -5.0
x_stop  = 5.0
x = np.arange(x_start, x_stop, dx)

# diabatic potentials
V_11 = V11(x, a, b)
V_22 = V22(x, a, b)

# adiabatic potentials for coupling c1=0.05 and c2=0.1
V_11_ad_1 = V11_adiabatic(x, a, b, c1, d)
V_22_ad_1 = V22_adiabatic(x, a, b, c1, d)

V_11_ad_2 = V11_adiabatic(x, a, b, c2, d)
V_22_ad_2 = V22_adiabatic(x, a, b, c2, d)

#calculate energy gap
delta_E_1 = min(V_11_ad_1 - V_22_ad_1)
delta_E_2 = min(V_11_ad_2 - V_22_ad_2)
print('Energy gap for c=0.05:\t ' + str(delta_E_1))
print('Energy gap for c=0.1:\t ' + str(delta_E_2))

# plot potentials
legend_fs = 16
axis_fs = 16
title_fs = 16

v11_color = 'C0'
v22_color = 'C1'

fig, ax = plt.subplots(1, 2, figsize=(15, 5))
ax[0].plot(x, V_11, color=v11_color, label=r'$V_{11}$')
ax[0].plot(x, V_22, color=v22_color, label=r'$V_{22}$')

ax[0].set_xlabel(r'x [Å]', fontsize=axis_fs)
ax[0].set_ylabel(r'Energy [eV]', fontsize=axis_fs)
ax[0].set_title(r'Diabatic potential energy', fontsize=title_fs)
ax[0].legend(fontsize=legend_fs)

ax[1].plot(x, V_11_ad_1, '--', color=v11_color, label=r'$V_{11}^{ad}$, $c = 0.05$')
ax[1].plot(x, V_11_ad_2, color=v11_color, label=r'$V_{11}^{ad}$, $c = 0.1$')
ax[1].plot(x, V_22_ad_1, '--', color=v22_color, label=r'$V_{22}^{ad}$, $c = 0.05$')
ax[1].plot(x, V_22_ad_2, color=v22_color, label=r'$V_{22}^{ad}$, $c = 0.1$')

ax[1].set_xlabel(r'x [Å]', fontsize=axis_fs)
ax[1].set_ylabel(r'Energy [eV]', fontsize=axis_fs)
ax[1].set_title(r'Adiabatic potential energy', fontsize=title_fs)
ax[1].legend(fontsize=legend_fs)

plt.show()
