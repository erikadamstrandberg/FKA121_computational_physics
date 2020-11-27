#%%
import numpy as np
import matplotlib.pyplot as plt

#%% 

Mb = int(100000)
B = int((N+1)/(Mb))

blocks = np.zeros([B,Mb])
f_blocks = np.zeros(B)

for i in range(len(blocks)):
    for j in range(len(f_blocks)):
        f_blocks[i] = f_blocks[i] + array[i+j]

print((1/B)*f_blocks)