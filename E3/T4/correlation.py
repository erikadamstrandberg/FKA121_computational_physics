#%%
import numpy as np
import matplotlib.pyplot as plt

#%% Auto correlation

array = np.genfromtxt('MC.txt', delimiter=',', skip_header=1)

N = len(array)

f  = ((1/N)*np.sum(array))**2
f2 =  (1/N)*np.sum(array**2)

#%% For k = 0 omega should be normalized

omega_k_0 = ((1/N)*np.dot(array,array) - f)/(f2 - f)
print(f"Should be normalized: {omega_k_0}")

forward = array - (1/N)*np.sum(array)
lag = array - (1/N)*np.sum(array)
omega_k_0 = (1/N)*np.dot(forward, lag)/(f2 - f)
print(f"Should be normalized: {omega_k_0}")

#%%
how_many = 200
ks = np.arange(0,how_many,1)
omega = np.zeros(how_many)
for k in ks:
    omega[k] = ((1/(N-k))*np.dot(array[0:N-k], array[k:N]) - f)/(f2 - f)

fig, ax = plt.subplots()
ax.plot(ks, omega)

#%%
how_many = 500
ks = np.arange(0,how_many,1)

omega = np.zeros(how_many)
for k in ks:
    forward = array[k:N] - (1/(N-k))*np.sum(array[k:N])
    lag = array[0:N-k] - (1/(N-k))*np.sum(array[0:N-k])
    omega[k] = (1/(N-k))*np.dot(lag, forward)/(f2 - f)
    
    
x = np.array([0, how_many])
y = np.array([0.13, 0.13])
    
fig, ax = plt.subplots()
ax.plot(ks, omega)
ax.plot(x, y)

ax.set_xlabel(r'$k$', fontsize='16')
ax.set_ylabel(r'$\Phi_k$', fontsize='16')
