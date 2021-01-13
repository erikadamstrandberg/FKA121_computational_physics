## Generating momentum for momentum operator
dp        = 2*np.pi*hbar_prim/(Nx*dx)              # Delta momentum
p_nyquist = 2*np.pi*hbar_prim/(2*dx)               # Nyquivst momentum
p         = np.arange(-p_nyquist, p_nyquist, dp)   # Generating momentum "freq"
p         = np.fft.fftshift(p)                     # Shifting momentum for FFT


V0    = 0.1                                        # Potential height [eV] 
alpha = 0.2                                        # Width [angstrom] 
V_x   = V0/(np.cosh(x/alpha)**2)                   # Eckart potential 

## Potential and momentum operator
v_prop =  np.exp(-1j*V_x*dt/hbar_prim)
p_prop = np.exp(-1j*p**2*dt/(hbar_prim*2*m_h))

## Function to propagate phi_x
def propagate(phi_x, p_prop, v_prop):
    return np.fft.ifft(p_prop*np.fft.fft(v_prop*phi_x))