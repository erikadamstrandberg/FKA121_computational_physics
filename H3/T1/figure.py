#%% H3b time-independent scroooodinger

# Standard includes
import numpy as np
import matplotlib.pyplot as plt

#%% Rescaling into ASU_prim units for H3b

u = 1.66054e-27    ## Atmoic weight in kg

T_asup = 1e-15      ## fs 
L_asup = 1e-10      ## angstrom
E_asup = 1.602e-19  ## eV

m_asup = E_asup*T_asup**2/L_asup**2 # Mass unit in ASU_prim
m_asup_u = m_asup/u                 # Mass unti in ASU_prim per u

hbar_si = 1.054571817e-34           # hbar in SI
hbar_asup = hbar_si*1e34/1.602      # hbar in ASU_prim. eV*fs

print(f"Mass in asup: {m_asup_u} u")
print(f"h_bar in asup: {hbar_asup} eV*fs")

#%% Plotting the Gaussian wave packet

def wave_package_position(x, x0, p0, d):
    hbar_si = 1.054571817e-34
    hbar_asup = hbar_si*1e34/1.602
    
    prefactor = np.sqrt(np.sqrt(np.pi)*d)
    envelope  = np.exp(-(x - x0)**2/(2*d**2))
    plane_wave = np.exp(1j*p0*(x - x0)/hbar_asup)
    return (1/prefactor)*envelope*plane_wave

def wave_package_position_evo(x, x0, t, p0, d, m):
    hbar_si = 1.054571817e-34
    hbar_asup = hbar_si*1e34/1.602
    
    prefactor = np.sqrt(np.sqrt(np.pi)*(d+1j*hbar_asup*t/(m*d)))
    envelope  = np.exp(-(x - x0 - p0*t/m)**2/(2*d**2*(1+1j*hbar_asup*t/(m*d**2))))
    plane_wave = np.exp(1j*p0*(x - x0 - p0*t/(2*m))/hbar_asup)
    return (1/prefactor)*envelope*plane_wave

d   = 0.5               # Width of our hydrogen atom
m_h = 1*m_asup_u        # Mass of our hydrogen atom

p0 = np.sqrt(0.2*m_h)   # Initial momentum of our hydrogen atom
x0 = 0.0                # Initial position of our hydrogen atom

dx      = 0.01
x_start = -10
x_stop  = 10
x = np.arange(x_start, x_stop, dx)
Nx = len(x)

phi_x = wave_package_position(x, x0, p0, d)
n_x = np.abs(phi_x)**2

phi_x_1 = wave_package_position_evo(x, x0, 0.0, p0, d, m_h)
n_x_1 = np.abs(phi_x_1)**2

fig, ax = plt.subplots()
ax.plot(x, n_x)
ax.plot(x, n_x_1)

plt.savefig('../T1/initial_gauss_position.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%% Fourier transforming to momentum space

dp = 1/(Nx*dx)                              # Delta momentum
p_nyquist = 1/(2*dx)                        # Nyquivst momentum
p = np.arange(-p_nyquist,p_nyquist,dp)      # Generating momentum "freq"

phi_p = np.fft.fftshift(np.fft.fft(phi_x))  # Shifting values from transform
n_p = np.abs(phi_p)**2                      # Probability density

x = np.array([0,0])
y = np.array([0,17500])

fig, ax = plt.subplots()
ax.plot(p, n_p)
ax.plot(x,y)

ax.set_xlim([-1,1])

plt.savefig('../T1/initial_gauss_momentum.pdf', format='pdf', bbox_inches='tight')
plt.show()


#%%

def create_V11(x, a, b):
    V11 = np.zeros(len(x))
    for i in range(len(x)):
        if (x[i] < 0):
            V11[i] = a*np.exp(x[i]/b)
        else:
            V11[i] = a*(2 - np.exp(-x[i]/b))
            
    return V11

def create_V22(x, a, b):
    return 2*a - create_V11(x, a, b)

def create_V12(x, c, d):
    return c*np.exp(-(x/d)**2)
 


x = np.linspace(-4,4,100)
a = 1
b = 1
c = 2
d = 1

V11 = create_V11(x, a, b)
V22 = create_V22(x, a, b)
V12 = create_V12(x, c, d)

fig, ax = plt.subplots()
ax.plot(x, V11)
ax.plot(x, V22)
ax.plot(x, V12)


#fig, ax = plt.subplots()
#ax.plot(x, np.real(phi_x))
#ax.plot(x, np.imag(phi_x)/np.max(np.imag(phi_x)))




#%% Uncorrelated! With only one electrion wavefunctions
signal = np.genfromtxt('../T1/signal.csv', delimiter=',')
time   = np.genfromtxt('../T1/time.csv', delimiter=',')
 

fig, ax = plt.subplots()
ax.plot(time, signal)


#%%

fft = np.fft.fft(signal)

#%%

spectrum   = np.genfromtxt('../T1/spectrum.csv', delimiter=',')
freq       = np.genfromtxt('../T1/freq.csv', delimiter=',')

fig, ax = plt.subplots()
ax.plot(freq, spectrum)

ax.set_xlim([0,5])

