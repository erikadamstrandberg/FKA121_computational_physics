#%%
import numpy as np
import matplotlib.pyplot as plt

#%%

array = np.genfromtxt('dist.csv', delimiter=',', skip_header=1)

uniform = array[1,:]