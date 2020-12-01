import numpy as np
import matplotlib.pyplot as plt

pos_data = np.genfromtxt('../T4/merged_q_trail.csv', delimiter=',')

x = pos_data[:, 0]
t = pos_data[:, -1]

plt.plot(t,x)
plt.show()

