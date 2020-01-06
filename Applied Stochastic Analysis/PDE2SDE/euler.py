import random
import numpy as np
from math import sqrt, pi, cos, sin
import matplotlib.pyplot as plt
from sklearn import linear_model
random.seed()
def euler(x, y, dt, N):
    out = i = 0
    while i < N:
        X = x
        Y = y
        while True:
            X += X * dt + random.gauss(0, sqrt(dt))
            Y += Y * dt + random.gauss(0, sqrt(dt))
            # notice that random.gauss(mu, sigma) means Normal(mu, sigma**2)
            if X**2 + Y**2 > 1:
                break
            else:
                out += (X**2 + Y**2 + 1) * dt
        i += 1
    return 1/2 - out / N


def main_a():
    Err = []; DT = []
    S_DT = [0.005*i for i in range(1, 31)]
    for sdt in S_DT:
        err = 0; dt = sdt**2; DT.append(dt)
        for i in range(20):
            th = 2 * pi * random.random()
            print(i)
            r = sqrt(random.random())
            x = r * cos(th); y = r * sin(th)
            err += euler(x, y, dt, 3000) - (x**2 + y**2)/2
        err /= 1
        Err.append(abs(err))
        print(err)

    print(Err)
    print(S_DT)

    S_DT = np.array(S_DT)
    Err = np.array(Err)
    log_DT = 2 * np.log(S_DT)
    log_Err = np.log(Err)
    regr = linear_model.LinearRegression()
    regr.fit(log_DT.reshape(-1, 1), log_Err)
    print('Score: %s' % regr.score(log_DT.reshape(-1, 1), log_Err))
    print('Coefficient: %s' % regr.coef_)
    plt.figure(figsize = (16, 9), dpi = 120)
    plt.plot(log_DT, log_Err, linewidth = 0, marker = 'o', markersize = 8)
    plt.plot(log_DT, regr.predict(log_DT.reshape(-1, 1)), color='red', linewidth=2)
    plt.ylabel("Log_Error")
    plt.xlabel("Log of Step Length")
    plt.show()

    regr = linear_model.LinearRegression()
    regr.fit(S_DT.reshape(-1, 1), Err)
    print('Score: %s' % regr.score(S_DT.reshape(-1, 1), Err))
    print('Coefficient: %s' % regr.coef_)
    plt.figure(figsize = (16, 9), dpi = 120)
    plt.plot(S_DT, Err, linewidth = 0, marker = 'o', markersize = 8)
    plt.plot(S_DT, regr.predict(S_DT.reshape(-1, 1)), color='red', linewidth=2)
    plt.ylabel("Error")
    plt.xlabel("Square Root of Step Length")
    plt.show()


main_a()
 