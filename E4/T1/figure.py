#%%
import numpy as np
import matplotlib.pyplot as plt

array = np.genfromtxt('../T1/gaussian.csv', delimiter=',')

step_bin = 0.1
start_bin = -5
stopp_bin = 5 + step_bin
b = np.arange(start_bin, stopp_bin, step_bin)

fig, ax = plt.subplots()
plt.hist(array, bins=b, density=True) 

