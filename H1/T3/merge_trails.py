import numpy as np
import sys

def main():
    nfiles = int(sys.argv[1])

    pos_data = np.genfromtxt('../T3/q_trail_0.csv', delimiter=',')
    vel_data = np.genfromtxt('../T3/v_trail_0.csv', delimiter=',')
    TPV = np.genfromtxt('../T3/TPV_0.csv', delimiter=',', skip_header=1)

    for i in range(nfiles-1):
        filenr = i+1

        pos_data_temp = np.genfromtxt('../T3/q_trail_'+ str(filenr) + '.csv', delimiter=',')
        pos_data = np.append(pos_data, pos_data_temp, axis=0)

        vel_data_temp = np.genfromtxt('../T3/v_trail_'+ str(filenr) + '.csv', delimiter=',')
        vel_data = np.append(vel_data, vel_data_temp, axis=0)

        TPV_temp = np.genfromtxt('../T3/TPV_'+ str(filenr) + '.csv', delimiter=',', skip_header=1)
        TPV = np.append(TPV, TPV_temp, axis=0)



    np.savetxt('merged_q_trail.csv', pos_data, delimiter=',') 
    np.savetxt('merged_v_trail.csv', vel_data, delimiter=',')
    np.savetxt('merged_TPV.csv', TPV, delimiter=',')

main()
