#%%
import numpy as np
import matplotlib.pyplot as plt

#%%
array = np.genfromtxt('../T2_avg/spectrum.csv', delimiter=',')
vfft = array[:,0]
freq = array[:,1]  

fig, ax = plt.subplots()
ax.plot(freq, vfft)

ax.set_xlim([0,10])


#%%
array = np.genfromtxt('../T2_avg/data/spectrum_high.csv', delimiter=',')
vfft_high = array[:,0]
freq = array[:,1]  

array = np.genfromtxt('../T2_avg/data/spectrum_low.csv', delimiter=',')
vfft_low = array[:,0]

fig, ax = plt.subplots()
ax.plot(freq, vfft_high, label=r'$P_{high} = 99.8kPa$')
ax.plot(freq, vfft_low,  label=r'$P_{low} = 2.75kPa$')

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

