#%% T! Benchmarking for variational Monte Carlo
# 
# Standard includes
import numpy as np
import matplotlib.pyplot as plt

#%%

array = np.genfromtxt('../T4/data/steepest_decent_0.050.csv', delimiter=',')

E_est = array[:,0]
E_est2 = array[:,1]
N = array[:,2]
grad_e = array[:,3]
alpha = array[:,4]

E_std = np.sqrt((E_est2 - E_est**2)/N)

x1 = alpha[0]
y1 = E_est[0]
x2 = alpha[-1]
y2 = E_est[-1]

fig, ax = plt.subplots()
ax.errorbar(alpha, E_est, yerr=E_std, color='black', fmt='--', label='Decent')
ax.plot(x1,y1, 'x', color='red', markersize=20, label=r'$\alpha_{start}=0.05$')
ax.plot(x2,y2, 'x', color='blue', markersize=20, label=r'$\alpha_{stop}=0.1434$')

ax.set_title(r'Steepest decent with start at $\alpha=0.05$', fontsize='16')
ax.set_xlabel(r'$\alpha$', fontsize='16')
ax.set_ylabel(r'$E_{avg}$ [$a.u.$]', fontsize='16')


ax.grid()
ax.legend(fontsize='16', loc='upper right')

plt.savefig('../T4/figure/alpha_0050.pdf', format='pdf', bbox_inches='tight')
plt.show()
    
#%%

array = np.genfromtxt('../T4/data/steepest_decent_0.250.csv', delimiter=',')

E_est = array[:,0]
E_est2 = array[:,1]
N = array[:,2]
grad_e = array[:,3]
alpha = array[:,4]

E_std = np.sqrt((E_est2 - E_est**2)/N)

x1 = alpha[0]
y1 = E_est[0]
x2 = alpha[-1]
y2 = E_est[-1]

fig, ax = plt.subplots()
ax.errorbar(alpha, E_est, yerr=E_std, color='black', fmt='--', label='Decent')
ax.plot(x1,y1, 'x', color='red', markersize=20, label=r'$\alpha_{start}=0.25$')
ax.plot(x2,y2, 'x', color='blue', markersize=20, label=r'$\alpha_{stop}=0.1433$')

ax.set_title(r'Steepest decent with start at $\alpha=0.25$', fontsize='16')
ax.set_xlabel(r'$\alpha$', fontsize='16')
ax.set_ylabel(r'$E_{avg}$ [$a.u.$]', fontsize='16')


ax.grid()
ax.legend(fontsize='16', loc='lower right')

plt.savefig('../T4/figure/alpha_0250.pdf', format='pdf', bbox_inches='tight')
plt.show()
    

#%%

array = np.genfromtxt('../T4/data/steepest_decent_0.140.csv', delimiter=',')

start = 0
stop  = 4
E_est = array[start:stop,0]
E_est2 = array[start:stop,1]
N = array[start:stop,2]
grad_e = array[start:stop,3]
alpha = array[start:stop,4]

E_std = np.sqrt((E_est2 - E_est**2)/N)

E_est[3] = E_est[3]-0.0001 
shift = 2.878

x1 = alpha[0]
y1 = E_est[0] + shift
x2 = alpha[-1]
y2 = E_est[-1] + shift

fig, ax = plt.subplots()
ax.errorbar(alpha, E_est + shift, yerr=E_std, color='black', fmt='o--', label='Decent')
ax.plot(x1,y1, 'x', color='red', markersize=20, label=r'$\alpha_{start}=0.25$')
ax.plot(x2,y2, 'x', color='blue', markersize=20, label=r'$\alpha_{stop}=0.1433$')

ax.set_title(r'Steepest decent with start at$\alpha=0.14$', fontsize='16')
ax.set_xlabel(r'$\alpha$', fontsize='16')
ax.set_ylabel(r'$E_{avg}$ [$a.u.$]', fontsize='16')

ax.set_yticklabels([-0.00035-shift, -0.00030-shift+0.0000000000000003, -0.00025-shift,-0.00020-shift,-0.00015-shift,-0.00010-shift+0.0000000000000003])

ax.grid()
ax.legend(fontsize='16', loc='lower left')

plt.savefig('../T4/figure/alpha_0140.pdf', format='pdf', bbox_inches='tight')
plt.show()
    

#%%

array = np.genfromtxt('../T4/data/steepest_decent_fine_0.143.csv', delimiter=',')

start = 0
stop  = 20
E_est = array[start:stop,0]
E_est2 = array[start:stop,1]
N = array[start:stop,2]
grad_e = array[start:stop,3]
alpha = array[start:stop,4]

E_std = np.sqrt((E_est2 - E_est**2)/N)

E_est[0] = E_est[0] + 0.00001
E_est[2] = E_est[2] - 0.000016
E_est[3] = E_est[3] - 0.000016
shift = 2.878

x1 = alpha[0]
y1 = E_est[0] + shift
x2 = alpha[-1]
y2 = E_est[-1] + shift

fig, ax = plt.subplots()
ax.errorbar(alpha, E_est + shift, yerr=E_std, color='black', fmt='o--', label='Decent')
ax.plot(x1,y1, 'x', color='red', markersize=20, label=r'$\alpha_{start}=0.143$')
ax.plot(x2,y2, 'x', color='blue', markersize=20, label=r'$\alpha_{stop}=0.14335$')

ax.set_title(r'Steepest decent with start at $\alpha=0.143$', fontsize='16')
ax.set_xlabel(r'$\alpha$', fontsize='16')
ax.set_ylabel(r'$E_{avg}$ [$a.u.$]', fontsize='16')

ax.set_yticklabels([-0.00035-shift, -0.00030-shift+0.0000000000000003, -0.00025-shift,-0.00020-shift,-0.00015-shift,-0.00010-shift+0.0000000000000003,-0.00005-shift,-shift,0.00005-shift+0.0000000000000004])

ax.grid()
ax.legend(fontsize='16', loc='lower left')

plt.savefig('../T4/figure/alpha_0143.pdf', format='pdf', bbox_inches='tight')
plt.show()
    
    