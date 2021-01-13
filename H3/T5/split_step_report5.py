## Generating momentum for momentum operator
dp        = 2*np.pi*hbar_prim/(Nx*dx)              # Delta momentum
p_nyquist = 2*np.pi*hbar_prim/(2*dx)               # Nyquivst momentum
p         = np.arange(-p_nyquist, p_nyquist, dp)   # Generating momentum "freq"
p         = np.fft.fftshift(p)                     # Shifting momentums for FFT 

## Creating momentum operator
P_prop = np.exp(-1j*p**2*dt/(hbar_prim*2*m_h))

## Generating adiabatic potentials 
V_a = V11_adiabatic(x, a, b, c, d, lim)     # V_a ground potental.
V_b = V22_adiabatic(x, a, b, c, d, lim)     # V_b excited potential.

## Creating potential operators
V_prop_1 = np.exp(-1j*V_a*dt/hbar_prim)
V_prop_2 = np.exp(-1j*V_b*dt/hbar_prim)

## Create transformation matrix
alpha = V_22 - V_11
beta = np.sqrt((V_11-V_22)**2 + 4*V_12**2)

A_11 = -(alpha+beta)/(2*V_12)
A_12 = -(alpha-beta)/(2*V_12)
A_21 = np.ones(len(x))
A_22 = np.ones(len(x))
A = np.array([[A_11, A_12], [A_21, A_22]])
A_inv = -(V_12/beta)*np.array([[A_22, -A_12], [-A_21, A_11]])

## Function to propagate phi_x
def propagate(phi_x, p_prop, v_prop_1, v_prop_2, A, A_inv, Nx):
    ## Initialize matrices
    phi_trans_v = np.zeros_like(phi_x)
    phi_trans_t = np.zeros_like(phi_x)
    phi_ad = np.zeros_like(phi_x)
    phi_dia = np.zeros_like(phi_x)
    
    ## Transformation to adiabatic representation to apply potential operator
    phi_ad = np.einsum('ijk, ik -> jk', A, phi_x)
       
    ## Applying potential operator
    phi_trans_v[0,:] = v_prop_1*phi_ad[0,:]
    phi_trans_v[1,:] = v_prop_2*phi_ad[1,:]
    
    ## Transformation to diabatic representation to apply kinetic operator
    phi_dia = np.einsum('ijk, ik -> jk', A_inv, phi_trans_v)
    
    ## Applying kinetic operator
    phi_trans_t[0,:] = np.fft.ifft(p_prop*np.fft.fft(phi_dia[0,:]))
    phi_trans_t[1,:] = np.fft.ifft(p_prop*np.fft.fft(phi_dia[1,:]))
    return phi_trans_t