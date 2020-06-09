# -*- coding: utf-8 -*-
"""
Created on Sat Jun  6 20:06:19 2020

@author: DELL
"""
import numpy as np
def Fehlberg_step(f, t, y, h):

    K1 = f(t, y)
    K2 = f(t + h / 4, y + h / 4 * K1)
    K3 = f(t + 3 * h / 8, y + (h * 3 / 32) * K1 + (h * 9 / 32) * K2)
    K4 = f(t + 12 * h / 13, y + (h * 1932 / 2197) * K1 - (h * 7200 / 2197) * K2 + (h * 7296 / 2197) * K3)
    K5 = f(t + h, y + (h * 439 / 216) * K1 - (h * 8) * K2 + (h * 3680 / 513) * K3 - (h * 845 / 4104) * K4)
    K6 = f(t + h / 2, y - (h * 8 / 27) * K1 + (h * 2) * K2 - (h * 3544 / 2565) * K3 + (h * 1859 / 4104) * K4 - (h * 11 / 40) * K5)
    
    y1 = y + h * (K1 * 16 / 135 + 6656 / 12825 * K3 + 28561 / 56430 * K4 - 9 / 50 * K5 + 2 / 55 * K6)
    y2 = y + h * (K1 * 25 / 216 + 1408 / 2565 * K3 + 2197 / 4104 * K4 - 1 / 5 * K5)

    return [y1, y2]


def Cash_Karp_step(f, t, y, h):
    
    K1 = f(t, y)
    K2 = f(t + h / 5, y + h / 5 * K1)
    K3 = f(t + 3 * h / 10, y + (h * 3 / 40) * K1 + (h * 9 / 40) * K2)
    K4 = f(t + 3 * h / 5, y + (h * 3 / 10) * K1 - (h * 9 / 10) * K2 + (h * 6 / 5) * K3)
    K5 = f(t + h, y - (h * 11 / 54) * K1 + (h * 5 / 2) * K2 - (h * 70 / 27) * K3 + (h * 35 / 27) * K4)
    K6 = f(t + 7 * h / 8, y + (h * 1631 / 55296) * K1 + (h * 175 / 512) * K2 + (h * 575 / 13824) * K3 + (h * 44275 / 110592) * K4 + (h * 253 / 4096) * K5)
    
    y1 = y + h * (K1 * 37 / 378 + 250 / 621 * K3 + 125 / 594 * K4 +  512 / 1771 * K6)
    y2 = y + h * (K1 * 2825 / 27648 + 18575 / 48384 * K3 + 13525 / 55296 * K4 + 277 / 14336 * K5 + 1 / 4 * K6)

    return [y1, y2]
        
        
def Dormand_Prince_step(f, t, y, h):
    
    K1 = f(t, y)
    K2 = f(t + h / 5, y + h / 5 * K1)
    K3 = f(t + 3 * h / 10, y + (h * 3 / 40) * K1 + (h * 9 / 40) * K2)
    K4 = f(t + 4 * h / 5, y + (h * 44 / 45) * K1 - (h * 56 / 15) * K2 + (h * 32 / 9) * K3)
    K5 = f(t + 8 / 9 * h, y + (h * 19372 / 6561) * K1 - (h * 25360 / 2187) * K2 + (h * 64448 / 6561) * K3 - (h * 212 / 729) * K4)
    K6 = f(t + h, y + (h * 9017 / 3168) * K1 - (h * 355 / 33) * K2 + (h * 46732 / 5247) * K3 + (h * 49 / 176) * K4 - (h * 5103 / 18656) * K5)
    K7 = f(t + h, y + (h * 35 / 384) * K1 + (h * 500 / 1113) * K3 + (h * 125 / 192) * K4 - (h * 2187 / 6784) * K5 + (h * 11 / 84) * K6)
    
    y1 = y + h * (K1 * 35 / 384 + 500 / 1113 * K3 + 125 / 192 * K4  -  2187 / 6784 * K5 + 11 / 84 * K6)
    y2 = y + h * (K1 * 5179 / 57600 + 7571 / 16695 * K3 + 393 / 640 * K4  -  92097 / 339200 * K5 + 187 / 2100 * K6 + 1 / 40 * K7)

    return [y1, y2]

def ode45_step(f, t, y, h, subfunc):
    
    SAFETY = 0.9
    PGROW = -0.2
    PSHRNK = -0.25
    ERRCON = (5 / SAFETY)**(1/PGROW)
    
    eps = 2**(-52)

    err_cal = np.abs(f(t, y)) + h * np.abs(f(t, y)) + eps
    out = subfunc(f, t, y, h)
    err = out[0]-out[1]
    errmax = max(np.abs(np.divide(err, err_cal)))/eps
    
    if errmax > 1:
        hnew = max(SAFETY * h * (errmax**PSHRNK), h / 10)
        tnew = t + hnew
        if tnew == t:
            hnew = h
            print('stepsize underflow')
    else:
        if errmax > ERRCON:
            hnew = SAFETY * h * (errmax**PGROW)
        else:
            hnew = 5 * h
    y = out[0]
    
    return [t + h, y, hnew]
    
    
def ode45(f, t, y, h, T, choice='Dormand_Prince'):
    if choice=='Dormand_Prince':
        subfunc = Dormand_Prince_step
    elif choice=='Cash_Karp':
        subfunc = Cash_Karp_step
    elif choice=='Fehlberg':
        subfunc = Fehlberg_step
    else:
        subfunc = Dormand_Prince_step
        print('incorrect subfunc, use Dormand_Prince instead')
    
    Y = [y]
    ts = [t]
    while t < T:
        [t, y, h] = ode45_step(f, t, y, h, subfunc)
        Y.append(y)
        ts.append(t)
    Y.pop(); ts.pop(); ts.append(T)
    t = ts[-2]
    h = T - t
    Y.append(Dormand_Prince_step(f, t, Y[-1], h)[0])
    return Y, ts

if __name__ == "__main__":

    def L(sigma, rho, beta):
        """
        Lorenz system
        """
        def f(t,x):
            f1 = sigma * (x[1] - x[0])
            f2 = x[0] * (rho - x[2]) - x[1]
            f3 = x[0] * x[1] - beta * x[2]
            return np.array([f1, f2, f3])
        return f
    
    
    sigma = 10; rho = 28; beta = 8 / 3
    f = L(sigma, rho, beta)
    h = 0.001
    T = 1
    y0 = np.array([1.0, 1.0, 1.0])
    Y0 = schemes_exp.Dormand_Prince(f, 0, y0, h, T)
    Y = ode54(f, 0, y0, h, T)
    res = np.linalg.norm(Y[0][-1]-Y0[-1])
    
    print(res)

