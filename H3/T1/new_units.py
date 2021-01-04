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
