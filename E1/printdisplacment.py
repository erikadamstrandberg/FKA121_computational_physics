#!/usr/bin/env python
###############################################################################
# E1code2
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv
#%%
array = np.genfromtxt('timetrail.csv', delimiter=',', skip_header=1)

# q1 = array[:, 0]
# q2 = array[:, 1]
# q3 = array[:, 2]
# U_kin = array[:, 3]
# U_pot = array[:, 4]
# time = array[:, 5]

fig, ax = plt.subplots()
ax.plot(array[:, 5], array[:, 0])
ax.plot(array[:, 5], array[:, 1])
ax.plot(array[:, 5], array[:, 2])

ax.set_xlabel('time (arb.unit)')
ax.set_ylabel('signal (arb.unit)')
ax.grid()

plt.show()
