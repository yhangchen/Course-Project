import numpy as np
from random import gauss

'''
mlmc(M, epsilon, DT, Richardson = True):
multi-level Monte Carlo path estimation
  M      = timestep refinement factor
  eps    = accuracy (rms error)
  mlmc_euler = function for level l estimator
output:
  value    = estimated value
  variance = variance of estimation
  L        = # of refinement layer

  mlmc_euler(layer, M, N, DT):
       l = level
       N = number of paths
       DT = initial time step
'''

def mlmc_euler(layer, M, N, DT):
    dt = DT / M**layer; dt_h = M * dt; sdt = np.sqrt(dt)
    out1 = out2 = out3 = out4 = 0; step = int(1/dt)
    for _ in range(N):
        X = 0; X_h = 0
        dWx_h = 0
        for count in range(1, step+1):
            dWx = gauss(0, sdt)
            X += - X * dt / 2 + dWx
            dWx_h +=  dWx
            if count % M == 0:
                X_h += - X_h * dt_h / 2 + dWx_h
                dWx_h = 0
        out1 += X**2 - X_h**2
        out2 += (X_h**2 - X**2)**2
        if not layer:
            out3 += X**2
            out4 += X**4
        
    if layer:
        return (0, out1, out2)
    else:
        return (0, out3, out4)




def mlmc(M, epsilon, DT, Richardson = True):
    L = -1; N = 100; stop = 0
    data0 = data1 = data2 = np.array([])
    while not stop:
        if L > 10:
            break
        L += 1
        data = mlmc_euler(L, M, N, DT)
        data0 = np.append(data0, N)
        data1 = np.append(data1, data[1])
        data2 = np.append(data2, data[2])
        print(data)
        V = np.divide(data2, data0) - (np.divide(data1, data0))**2
        M_L = np.power(M, np.array(range(L+1)))
        N_L = np.ceil(2/(epsilon**2) * np.sqrt(np.divide(V, M_L)) * np.sum(np.sqrt(np.multiply(V, M_L))))
        diff_N_L = N_L - data0
        for i in range(L+1):
            if diff_N_L[i] > 0:
                data = mlmc_euler(i, M, int(diff_N_L[i]), DT)
                data0[i] += diff_N_L[i]
                data1[i] += data[1]
                data2[i] += data[2]

        Y_L = data2[L]/data1[L]; Y_L_1 = 1/M * data2[L-1]/data1[L-1]
        if Richardson:
            if L > 1 and M**L > 50:
                stop = int(abs(Y_L - Y_L_1) < (M**2-1)*epsilon/np.sqrt(2)) + int(M**L > 100000)
        else:
            if L > 1 and M**L > 50:
                stop = int(max(abs(Y_L), abs(Y_L_1)) < (M-1)*epsilon/np.sqrt(2)) + int(M**L > 100000)
    value = np.sum(np.divide(data1, data0))
    variance = np.sum(np.divide(data2, data0) - (np.divide(data1, data0))**2)
    if Richardson:
        value += data1[L]/data0[L]/(M-1)
    return value, variance, L


a = mlmc(3, 0.005, 0.01, Richardson = True)
print(a)
print(a[0]-1+1/np.e)