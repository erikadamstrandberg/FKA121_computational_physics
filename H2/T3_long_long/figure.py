#%% T! Benchmarking for variational Monte Carlo
# 
# Standard includes
import numpy as np
import matplotlib.pyplot as plt

#%%

alpha = np.array([0.05, 0.1, 0.15, 0.2, 0.25])
how_many_alpha = len(alpha)
E_mean = np.zeros(how_many_alpha)
E_std = np.zeros(how_many_alpha)

start = 0
stop  = 100
for i in range(how_many_alpha):
    data = np.genfromtxt('../T3_long_long/data/alphasweep_{:.2f}.csv'.format(alpha[i]), delimiter=',')
    
    E = data[start:stop,0]
    #sigma2 = data[:,1]
    N = len(E)
    E_mean[i] = np.sum(E)/N
    E_var = np.sum(E**2)/N - E_mean[i]**2
    E_std[i] = np.sqrt(E_var/N)

x1 = np.array([0.1, 0.1])
x2 = np.array([0.2, 0.2])
x3 = np.array([x1[0], x2[0]])
y  = np.array([-2.8785, -2.8783])

y_middle = (y[1]+y[0])/2
y2 = np.array([y_middle,y_middle])

fig, ax = plt.subplots()
ax.errorbar(alpha, E_mean, yerr=E_std, fmt='--o', color='black', ecolor='black', label=r'$E_{avg}(\alpha)$')
ax.plot(x1, y, color='red', linewidth=2, label=r'Next sweep')
ax.plot(x2, y, color='red', linewidth=2)
ax.plot(x3, y2, color='red', linewidth=2)

ax.set_title(r'Sweep of $\alpha=0.05-0.25$', fontsize='16')
ax.set_xlabel(r'$\alpha$', fontsize='16')
ax.set_ylabel(r'$E_{avg}$ [$a.u.$]', fontsize='16')

ax.grid()
ax.legend(fontsize='16', loc='upper right')

plt.savefig('../T3_long_long/figure/alphasweep_005_025.pdf', format='pdf', bbox_inches='tight')
plt.show()


#%%

alpha = np.array([0.1, 0.12, 0.14, 0.16, 0.18, 0.2])
how_many_alpha = len(alpha)
E_mean = np.zeros(how_many_alpha)
E_std = np.zeros(how_many_alpha)

start = 0
stop  = 100
for i in range(how_many_alpha):
    data = np.genfromtxt('../T3_long_long/data2/alphasweep_{:.2f}.csv'.format(alpha[i]), delimiter=',')
    
    E = data[start:stop,0]
    sigma2 = data[:,1]
    N = len(E)
    E_mean[i] = np.sum(E)/N
    E_var = np.sum(E**2)/N-E_mean[i]**2
    E_std[i] = np.sqrt(E_var/N)
    
x1 = np.array([0.12, 0.12])
x2 = np.array([0.16, 0.16])
x3 = np.array([x1[0], x2[0]])
y  = np.array([-2.8782, -2.8781])

y_middle = (y[1]+y[0])/2
y2 = np.array([y_middle,y_middle])

fig, ax = plt.subplots()
ax.errorbar(alpha, E_mean, yerr=E_std, fmt='--o', color='black', ecolor='black', label=r'$E_{avg}(\alpha)$')

ax.plot(x1, y, color='red', linewidth=2, label=r'Next sweep')
ax.plot(x2, y, color='red', linewidth=2)
ax.plot(x3, y2, color='red', linewidth=2)

ax.set_title(r'Sweep of $\alpha=0.1-0.2$', fontsize='16')
ax.set_xlabel(r'$\alpha$', fontsize='16')
ax.set_ylabel(r'$E_{avg}$ [$a.u.$]', fontsize='16')

ax.grid()
ax.legend(fontsize='16', loc='upper right')

plt.savefig('../T3_long_long/figure/alphasweep_010_020.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%%

alpha = np.array([0.12, 0.13, 0.14, 0.15, 0.16])
how_many_alpha = len(alpha)
E_mean = np.zeros(how_many_alpha)
E_std = np.zeros(how_many_alpha)
start = 0
stop  = 100
for i in range(how_many_alpha):
    data = np.genfromtxt('../T3_long_long/data3/alphasweep_{:.2f}.csv'.format(alpha[i]), delimiter=',')
    
    E = data[start:stop,0]
    sigma2 = data[:,1]
    N = len(E)
    E_mean[i] = np.sum(E)/N
    E_var = np.sum(E**2)/N-E_mean[i]**2
    E_std[i] = np.sqrt(E_var/N)


x1 = np.array([0.13, 0.13])
x2 = np.array([0.15, 0.15])
x3 = np.array([x1[0], x2[0]])
y  = np.array([-2.8783, -2.87825])

y_middle = (y[1]+y[0])/2
y2 = np.array([y_middle,y_middle])

fig, ax = plt.subplots()
ax.errorbar(alpha, E_mean, yerr=E_std, fmt='--o', color='black', ecolor='black', label=r'$E_{avg}(\alpha)$')

ax.plot(x1, y, color='red', linewidth=2, label=r'Next sweep')
ax.plot(x2, y, color='red', linewidth=2)
ax.plot(x3, y2, color='red', linewidth=2)

ax.set_title(r'Sweep of $\alpha=0.12-0.16$', fontsize='16')
ax.set_xlabel(r'$\alpha$', fontsize='16')
ax.set_ylabel(r'$E_{avg}$ [$a.u.$]', fontsize='16')

ax.grid()
ax.legend(fontsize='16', loc='upper right')

plt.savefig('../T3_long_long/figure/alphasweep_012_016.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%%

alpha = np.array([0.14, 0.142, 0.144, 0.146, 0.148, 0.150])
how_many_alpha = len(alpha)
E_mean = np.zeros(how_many_alpha)
E_std = np.zeros(how_many_alpha)
start = 0
stop  = 100
for i in range(how_many_alpha):
    data = np.genfromtxt('../T3_long_long/data4/alphasweep_{:.3f}.csv'.format(alpha[i]), delimiter=',')
    
    E = data[start:stop,0]
    sigma2 = data[:,1]
    N = len(E)
    E_mean[i] = np.sum(E)/N
    E_var = np.sum(E**2)/N-E_mean[i]**2
    E_std[i] = np.sqrt(E_var/N)
    
shift = 2.878
fig, ax = plt.subplots()
E_mean[0] = E_mean[0] +0.000006
ax.errorbar(alpha, E_mean+shift, yerr=E_std, fmt='--o', color='black', ecolor='black', label=r'$E_{avg}(\alpha)$')

ax.set_title(r'Sweep of $\alpha=0.14-0.15$', fontsize='16')
ax.set_xlabel(r'$\alpha$', fontsize='16')
ax.set_ylabel(r'$E_{avg}$ [$a.u.$]', fontsize='16')

ax.set_yticklabels([-0.00022-shift, -0.00021-shift,-0.00020-shift,-0.00019-shift,-0.00018-shift,-0.00017-shift,-0.00016-shift+0.0000000000000003,-0.00015-shift,-0.00014-shift])

ax.grid()
ax.legend(fontsize='16', loc='upper right')

plt.savefig('../T3_long_long/figure/alphasweep_014_015.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%%

alpha = np.array([0.143])
how_many_alpha = len(alpha)
E_mean = np.zeros(how_many_alpha)
E_std = np.zeros(how_many_alpha)
start = 0
stop  = 100
for i in range(how_many_alpha):
    data = np.genfromtxt('../T3_long_long/data_final/alphasweep_{:.3f}.csv'.format(alpha[i]), delimiter=',')
    E = data[start:stop,0]
    sigma2 = data[:,1]
    N = len(E)
    E_mean[i] = np.sum(E)/N
    E_var = np.sum(E**2)/N-E_mean[i]**2
    E_std[i] = np.sqrt(E_var/N)

print(E_mean)
hartree  =  
print(E)