## Generating momentum for momentum operator
dp        = 2*np.pi*hbar_prim/(Nx*dx)              # Delta momentum
p_nyquist = 2*np.pi*hbar_prim/(2*dx)               # Nyquivst momentum
p         = np.arange(-p_nyquist, p_nyquist, dp)   # Generating momentum "freq"
p         = np.fft.fftshift(p)                     # Shifting to match FFT

## Operator for free space propagation
V_prop = 1
P_prop = np.exp(-1j*p**2*dt/(hbar_prim*2*m_h))

## Function to propagate phi_x in free space
def propagate(phi_x, p_prop, v_prop):
    return np.fft.ifft(p_prop*np.fft.fft(v_prop*phi_x))