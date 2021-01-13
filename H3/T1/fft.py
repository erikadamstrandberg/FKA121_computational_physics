## Generating the Gaussian wave packet
def wave_package_position(x, x0, p0, d, hbar_prim):    
    prefactor = 1/np.sqrt(np.sqrt(np.pi)*d)
    envelope  = np.exp(-(x - x0)**2/(2*d**2))
    plane_wave = np.exp(1j*p0*(x - x0)/hbar_prim)
    return prefactor*envelope*plane_wave

## Gaussian wave packet in momentum space
def wave_package_momentum(p, p0, d, hbar_prim):
    prefactor = np.sqrt(d/(np.sqrt(np.pi)*hbar_prim))
    envelope  = np.exp(-d**2*(p0 - p)**2/(2*hbar_prim**2))
    return prefactor*envelope


phi_x = wave_package_position(x, x0, p0, d, hbar_prim)  # Generate wave packet
n_x   = np.abs(phi_x)**2                                # Proability density

dp        = 2*np.pi*hbar_prim/(Nx*dx)                   # Delta momentum
p_nyquist = 2*np.pi*hbar_prim/(2*dx)                    # Nyquivst momentum
p         = np.arange(-p_nyquist, p_nyquist, dp)        # Generating momentum "freq"

phi_p = np.fft.fftshift(np.fft.fft(phi_x))              # Shifting values from transform
n_p   = np.abs(phi_p)**2*dx**2/(2*np.pi*hbar_prim)      # Probability density with norm