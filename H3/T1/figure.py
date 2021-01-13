#%% H3b time-independent scroooodinger

# Standard includes
import numpy as np
import matplotlib.pyplot as plt

#%% Rescaling into ASU_prim units for H3b

u = 1.66054e-27    ## Atmoic weight in kg

T_prim = 1e-15      ## fs 
L_prim = 1e-10      ## angstrom
E_prim = 1.602e-19  ## eV

m_prim = E_prim*T_prim**2/L_prim**2 # Mass unit in ASU_prim
m_prim_u = m_prim/u                 # Mass unti in ASU_prim per u

hbar_si   = 1.054571817e-34         # hbar in SI
hbar_prim = hbar_si/(E_prim*T_prim) # hbar in ASU_prim. eV*fs

print(f"Mass in asup: {m_prim_u} u")
print(f"h_bar in asup: {hbar_prim} eV*fs")

#%% Plotting the Gaussian wave packet

def wave_package_position(x, x0, p0, d, hbar_prim):    
    prefactor = np.sqrt(np.sqrt(np.pi)*d)
    envelope  = np.exp(-(x - x0)**2/(2*d**2))
    plane_wave = np.exp(1j*p0*(x - x0)/hbar_prim)
    return (1/prefactor)*envelope*plane_wave

def wave_package_position_evo(x, x0, t, p0, d, m, hbar_prim):   
    prefactor = np.sqrt(np.sqrt(np.pi)*(d+1j*hbar_prim*t/(m*d)))
    envelope  = np.exp(-(x - x0 - p0*t/m)**2/(2*d**2*(1 + 1j*hbar_prim*t/(m*d**2))))
    plane_wave = np.exp(1j*p0*(x - x0 - p0*t/(2*m))/hbar_prim)
    return (1/prefactor)*envelope*plane_wave

def x_width_anal(t, d, m, hbar_prim):
    prefactor = d/np.sqrt(2)
    time_factor = np.sqrt(1 + hbar_prim**2*t**2/(m**2*d**4))
    return prefactor*time_factor

def p_width_anal(t, d, m, hbar_prim):
    prefactor = d/np.sqrt(2)
    time_factor = np.sqrt(1 + hbar_prim**2*t**2/(m**2*d**4))
    return prefactor*time_factor

def std_from_FWHM(n, dx):
    max_value = np.max(n)
    width_limit = max_value/2
    return np.sum(n > width_limit)*dx/(2*np.sqrt(2*np.log(2)))

d   = 0.5               # Width of our hydrogen atom
m_h = 1/m_prim_u        # Mass of our hydrogen atom

p0 = np.sqrt(0.2*m_h)  # Initial momentum of our hydrogen atom
x0 = 3.0               # Initial position of our hydrogen atom

dx      = 0.0001
x_start = -50
x_stop  = 50.0
x       = np.arange(x_start, x_stop, dx)
Nx      = len(x)

phi_x = wave_package_position(x, x0, p0, d, hbar_prim)
n_x   = np.abs(phi_x)**2

t1      = 20
phi_x_1 = wave_package_position_evo(x, x0, t1, p0, d, m_h, hbar_prim)
n_x_1   = np.abs(phi_x_1)**2

width = std_from_FWHM(n_x_1, dx)
width_anal = x_width_anal(t1, d, m_h, hbar_prim)

total_probability = np.sum(n_x_1*dx)

print(f'Width from finding the FWHM of the wavepacket: {width} Å')
print(f'Analytical width:                              {width_anal} Å')
print(f'Total probability: {total_probability}')

plot_x0_x = np.array([x0,x0])
y_min = -0.01
y_max = np.max(n_x)+0.1
plot_x0_y = np.array([y_min,y_max])

fig, ax = plt.subplots()
ax.plot(x, n_x, color='blue', linewidth=4, label=r'$|\psi(x)|^2$')
ax.plot(plot_x0_x, plot_x0_y, color='black', label=r'$x_0=3.0\ Å$')

ax.set_title(r'Position of initial Gaussian wave packet', fontsize='16')
ax.set_xlabel(r'$x$ [$Å$]', fontsize='16')
ax.set_ylabel(r'$|\psi(x)|^2$ [$Å^{-1}$]', fontsize='16')

ax.grid()
ax.legend(fontsize='16', loc='upper right')

ax.set_xlim([0,6])
ax.set_ylim([y_min,y_max])

plt.savefig('../T1/initial_gauss_position.pdf', format='pdf', bbox_inches='tight')
plt.show()

#%% Fourier transforming to momentum space

def wave_package_momentum(p, p0, d, hbar_asup):
    prefactor = np.sqrt(d/(np.sqrt(np.pi)*hbar_asup))
    envelope  = np.exp(-d**2*(p0 - p)**2/(2*hbar_asup**2))
    return prefactor*envelope

dp        = 2*np.pi*hbar_prim/(Nx*dx)              # Delta momentum
p_nyquist = 2*np.pi*hbar_prim/(2*dx)               # Nyquivst momentum
p         = np.arange(-p_nyquist, p_nyquist, dp)   # Generating momentum "freq"

phi_p = np.fft.fftshift(np.fft.fft(phi_x))         # Shifting values from transform
n_p   = np.abs(phi_p)**2*dx**2/(2*np.pi*hbar_prim) # Probability density

print(np.sum(n_p)*dp)

phi_p_anal = wave_package_momentum(p, p0, d, hbar_prim)
n_p_anal = np.abs(phi_p_anal)**2

print(np.sum(n_p_anal)*dp)

print(std_from_FWHM(n_p, dp))
print(hbar_prim/(np.sqrt(2)*d))

k = p/hbar_prim
k0 = p0/hbar_prim

plot_p0_x = np.array([k0, k0])
y_min = -0.01
y_max = np.max(n_p)+0.1
plot_p0_y = np.array([y_min, y_max])

fig, ax = plt.subplots()
ax.plot(k, n_p, color='lightgreen', linewidth=4      , label=r'$|\psi(k)|^2$')
ax.plot(k, n_p_anal, '--', color='black', linewidth=4, label=r'$|\psi(k)|^2_{analytical}$')
ax.plot(plot_p0_x, plot_p0_y, color='black'          , label=r'$k_0=6.91\ Å^{-1}$')

ax.set_title(r'Momentum of initial Gaussian wave packet', fontsize='16')
ax.set_xlabel(r'$k$ [$Å^{-1}$]', fontsize='16')
ax.set_ylabel(r'$|\psi(k)|^2$ [$Å$]', fontsize='16')

ax.grid()
ax.legend(fontsize='16', loc='upper right')

ax.set_xlim([0,14])
ax.set_ylim([y_min,y_max])

#plt.savefig('../T1/initial_gauss_momentum.pdf', format='pdf', bbox_inches='tight')
#plt.show()


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
#ax.plot(x, V12)


#fig, ax = plt.subplots()
#ax.plot(x, np.real(phi_x))
#ax.plot(x, np.imag(phi_x)/np.max(np.imag(phi_x)))



