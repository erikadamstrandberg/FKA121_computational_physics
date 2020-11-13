import numpy as np

def main():
    # SI units
    m0 = 12 # kg
    v0 = 3  # m/s
    k0 = 1000 # kg/s^2
    l0 = 4  # m

    E0_kin = m0*v0**2/2
    E0_pot = k0*l0**2/2

    # what we want to have in our units
    mp = 1  # --> m = 1/m0
    kp = 1  # --> k = 1/k0
    
    # scalefactors
    m = 1/m0
    k = 1/k0
    t = np.sqrt(m/k)
    l = 1e10        # (can choose this arbitrarily)

    E = m*l**2/t**2
    v = l/t
    
    # transformed quantities 
    vp = v*v0
    lp = l*l0
    # Energy calculated in our units
    Ep_kin = mp*vp**2/2
    Ep_pot = kp*lp**2/2

    # transform back to SI units 
    E0_kin_check = Ep_kin/E
    E0_pot_check = Ep_pot/E


    print(E0_pot)
    print(E0_pot_check)

main()


    


