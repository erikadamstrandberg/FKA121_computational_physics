#%%
import numpy as np
import matplotlib.pyplot as plt

#%% Auto correlation

array = np.genfromtxt('MC.txt', delimiter=',', skip_header=1)

N = len(array)

f  = ((1/N)*sum(array))**2
f2 =  (1/N)*sum(array**2)


#%% For k = 0 omega should be normalized

omega_k_0 = ((1/N)*sum(array*array) - f)/(f2 - f)
print(omega_k_0)

#%% Sweep of omegea!

how_many = 3
ks = np.arange(0,how_many,1)

omega = np.zeros(how_many)

for k in ks:
    omega[k] = ((1/(N-k))*sum(array[0:N-k]*array[k:N]) - f)/(f2 - f)
    
fig, ax = plt.subplots()
ax.plot(ks, omega)

#%% 

Mb = int(100000)
B = int((N+1)/(Mb))

blocks = np.zeros([B,Mb])
f_blocks = np.zeros(B)

for i in range(len(blocks)):
    for j in range(len(f_blocks)):
        f_blocks[i] = f_blocks[i] + array[i+j]

print((1/B)*f_blocks)