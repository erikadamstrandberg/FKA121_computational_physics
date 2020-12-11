#%%
import numpy as np
import scipy
import matplotlib.pyplot as plt

#%%
array = np.genfromtxt('../T4/correlation_low.csv', delimiter=',')
corr = array[:,0]
time = array[:,1] 
dt = time[1] - time[0] 
N = len(corr)

corr_inv = np.zeros(N)
N2 = int(N/2)
corr_half = corr[0:N2]

phi_fft = np.fft.fft(corr_half)
freq = np.fft.fftfreq(N/2, d=dt)
fig, ax = plt.subplots()
ax.plot(freq, phi_fft) 

#ax.set_xlim([0,10])






#%%
array = np.genfromtxt('../T4/data/power_corr_low.csv', delimiter=',')
spectrum_low = array[:,0]
freq = array[:,1]

array = np.genfromtxt('../T4/data/power_corr_high.csv', delimiter=',')
spectrum_high = array[:,0]

fig, ax = plt.subplots()
ax.plot(freq, spectrum_high,label=r'$P_{high} = 99.8kPa$')
ax.plot(freq, spectrum_low, label=r'$P_{low} = 2.75kPa$')

ax.set_xlim([0,10])

ax.set_title(r'Spectrum from correlation function', fontsize='16')
ax.set_xlabel(r'$f$ [$kHz$]', fontsize='16')
ax.set_ylabel(r'$P$ [$arb.u.$]', fontsize='16')

ax.grid()
ax.legend(fontsize='16', loc='upper right')
plt.show()