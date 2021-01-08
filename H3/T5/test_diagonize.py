import numpy as np
import matplotlib.pyplot as plt

#%%
M = np.array([[3, 2], [1, 4]])

a = np.array([1.0,1.0,1.0,1.0,1.0,1.0])
b = np.array([1.0,1.0,1.0,1.0,1.0,1.0])
c = np.array([-2.0,-2.0,-2.0,-2.0,-2.0,-2.0])
d = np.array([1.0,1.0,1.0,1.0,1.0,1.0])
A = np.array([[c, b], [a, d]])

#A_inv = np.array([[-1, 1],[1, 2]])/3
a = np.array([-1.0,-1.0,-1.0,-1.0,-1.0,-1.0])/3.0
b = np.array([1.0,1.0,1.0,1.0,1.0,1.0])/3.0
c = np.array([1.0,1.0,1.0,1.0,1.0,1.0])/3.0
d = np.array([2.0,2.0,2.0,2.0,2.0,2.0])/3.0
A_inv = np.array([[a, b], [c, d]])

x = np.array([10.0, 4.0, 2.0, 2.0, 5.0, 11.0])
y = np.array([-4.0,-100.0,-11.0, -2.0, -10.0, 20.0])
phi = np.array([x, y])
phi_trans = np.zeros_like(phi)
phi_back = np.zeros_like(phi)

    
for i in range(6):
    print(i)
    phi_trans[:,i] = A[:,:,i] @ phi[:,i]
    phi_back[:,i]  = A[:,:,i] @ A_inv[:,:,i] @ phi[:,i]
    


einsum_trans = np.einsum('ijk,ik->jk', A, phi)
einsum_back  = np.einsum('ijk,ik->jk', A_inv, einsum_trans)
    
print(phi)
print(einsum_back)  
