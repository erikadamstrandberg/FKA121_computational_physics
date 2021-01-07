import numpy as np
import matplotlib.pyplot as plt

M = np.array([[3, 2], [1, 4]])

v1 = np.array([1, 1])
v2 = np.array([-2, 1])
#A = np.array([v2, v1]).T
a = np.array([1,1,1,1,1])
b = np.array([1,1,1,1,1])
c = np.array([-2,-2,-2,-2,-2])
d = np.array([1,1,1,1,1])
A = np.array([[c, b], [a, d]])

#A_inv = np.array([[-1, 1],[1, 2]])/3
a = np.array([-1,-1,-1,-1,-1])/3
b = np.array([1,1,1,1,1])/3
c = np.array([1,1,1,1,1])/3
d = np.array([2,2,2,2,2])/3
A_inv = np.array([[a, b], [c, d]])


#trans = A_inv @ M @ A
#back = A @ trans @ A_inv
#print(back)



x = np.array([2, 3, 4, 5, 6])
y = np.array([-4,-3,-2, -1, -5])
phi = np.array([x, y])
phi_trans = np.zeros_like(phi)
phi_back = np.zeros_like(phi).T



for i in range(5):
    print(i)
    A_temp = A[:,:,i]
    A_inv_temp = A_inv[:,:,i]
    phi_temp = phi[:,i]

    phi_trans_temp = A_temp @ phi_temp
    phi_trans[:,i] = phi_trans_temp
    #print(phi_trans_temp)
    #print(phi_trans)
    #print(A_inv_temp @ phi_trans_temp)
    print(phi[:,i])
    print(A_inv_temp @ phi_trans[:,i])
    phi_back_temp =  A_inv_temp @ phi_trans[:,i]
    print(phi_back_temp)
    phi_back[i,:] = phi_back_temp.T
    print(phi_back[i,:])

    #phi_trans[:,i] = A[:,:,i] @ phi[:,i]
    #print(phi_back_temp)
    #phi_back[:,i] = A_inv[:,:,i] @ phi_trans[:,i]
    #print(phi_back[:,i])

#print(phi_trans)
#print(phi)
#print(phi_back)

