#%% T! Benchmarking for variational Monte Carlo
# 
# Standard includes
import numpy as np
import matplotlib.pyplot as plt

#%%

array = np.genfromtxt('../T3/data/finalvalue_005.csv', delimiter=',')

E_l = array[0]
sigma2 = array[1]
alpha = array[2]
N = array[3]
ns = array[4]

pm = np.sqrt(sigma2)/(np.sqrt(N/ns))

print(f"E_l = <E_l> +- sigma/sqrt(N/ns) for alpha = {alpha}")
print(f"E_l = {E_l} +- {pm}")

#%%

array = np.genfromtxt('../T3/data/finalvalue_010.csv', delimiter=',')

E_l = array[0]
sigma2 = array[1]
alpha = array[2]
N = array[3]
ns = array[4]

pm = np.sqrt(sigma2)/(np.sqrt(N/ns))

print(f"E_l = <E_l> +- sigma/sqrt(N/ns) for alpha = {alpha}")
print(f"E_l = {E_l} +- {pm}")

#%%

array = np.genfromtxt('../T3/data/finalvalue_015.csv', delimiter=',')

E_l = array[0]
sigma2 = array[1]
alpha = array[2]
N = array[3]
ns = array[4]

pm = np.sqrt(sigma2)/(np.sqrt(N/ns))

print(f"E_l = <E_l> +- sigma/sqrt(N/ns) for alpha = {alpha}")
print(f"E_l = {E_l} +- {pm}")

#%%

array = np.genfromtxt('../T3/data/finalvalue_020.csv', delimiter=',')

E_l = array[0]
sigma2 = array[1]
alpha = array[2]
N = array[3]
ns = array[4]

pm = np.sqrt(sigma2)/(np.sqrt(N/ns))

print(f"E_l = <E_l> +- sigma/sqrt(N/ns) for alpha = {alpha}")
print(f"E_l = {E_l} +- {pm}")

#%%

array = np.genfromtxt('../T3/data/finalvalue_025.csv', delimiter=',')

E_l = array[0]
sigma2 = array[1]
alpha = array[2]
N = array[3]
ns = array[4]

pm = np.sqrt(sigma2)/(np.sqrt(N/ns))

print(f"E_l = <E_l> +- sigma/sqrt(N/ns) for alpha = {alpha}")
print(f"E_l = {E_l} +- {pm}")
