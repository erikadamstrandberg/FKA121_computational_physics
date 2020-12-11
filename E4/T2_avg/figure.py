#%%
import numpy as np
import matplotlib.pyplot as plt

#%%

kb = 1.38e-23
T = 297
rho = 2650
d = 2.79e-6
r = d/2
m_pearl = rho*4*np.pi*r**3/3
w0 = 3.1*2*np.pi

tau_low = 147.3*1e-3
tau_high = 48.5*1e-3

mu_low = 1/tau_low
mu_high = 1/tau_high

f = np.linspace(0,20,1000)
omega = 2*np.pi*f
anal_low = (kb*T/m_pearl)*(2*mu_low*omega**2/((omega**2-w0**2)**2+mu_low**2*omega**2))
anal_high = (kb*T/m_pearl)*(2*mu_high*omega**2/((omega**2-w0**2)**2+mu_high**2*omega**2))


array = np.genfromtxt('../T2_avg/data/spectrum_high.csv', delimiter=',')
vfft_high = array[:,0]
freq = array[:,1]  

array = np.genfromtxt('../T2_avg/data/spectrum_low.csv', delimiter=',')
vfft_low = array[:,0]

fig, ax = plt.subplots()
ax.plot(freq, vfft_high/0.001, label=r'$P_{high} = 99.8kPa$')
ax.plot(freq, vfft_low/0.001,  label=r'$P_{low} = 2.75kPa$')
#ax.plot(f, anal_low*1.3e4,  '--', label=r'$P_{a,low}$', linewidth=3)
#ax.plot(f, anal_high*1.4e4, '--', label=r'$P_{a,high}$', linewidth=3)

ax.set_xlim([0,10])

ax.set_title(r'Averaged spectrum from v-trail, $d\tau = 0.05\ ms$', fontsize='16')
ax.set_xlabel(r'$f$ [$kHz$]', fontsize='16')
ax.set_ylabel(r'$P$ [$arb.u.$]', fontsize='16')


ax.grid()
ax.legend(fontsize='16', loc='upper right')
plt.show()


#%%
array = np.genfromtxt('../T2_avg/data/spectrum_high_fine.csv', delimiter=',')
vfft_high = array[:,0]
freq = array[:,1]  

array = np.genfromtxt('../T2_avg/data/spectrum_low_fine.csv', delimiter=',')
vfft_low = array[:,0]

fig, ax = plt.subplots()
ax.plot(freq, vfft_high, label=r'$P_{high} = 99.8kPa$')
ax.plot(freq, vfft_low,  label=r'$P_{low} = 2.75kPa$')

ax.set_xlim([0,10])

ax.set_title(r'Averaged spectrum from v-trail, $d\tau = 0.025\ ms$', fontsize='16')
ax.set_xlabel(r'$f$ [$kHz$]', fontsize='16')
ax.set_ylabel(r'$P$ [$arb.u.$]', fontsize='16')

ax.grid()
ax.legend(fontsize='16', loc='upper right')
plt.show()

