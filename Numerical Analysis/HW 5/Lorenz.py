# /*
#  * @Author: Yihang Chen 
#  * @Date: 2020-05-01 10:31:21 
#  * @Last Modified by: Yihang Chen
#  * @Last Modified time: 2020-05-15 10:40:32
#  */

import numpy as np, matplotlib.pyplot as plt, glob, os
import IPython.display as IPdisplay, matplotlib.font_manager as fm
from mpl_toolkits.mplot3d.axes3d import Axes3D
from PIL import Image


def RungeKutta(f, t, y, h, T):
    Y = [y]

    while t < T:
        K1 = f(t, y)
        K2 = f(t + h / 2, y + (h / 2) * K1)
        K3 = f(t + h / 2, y + (h / 2) * K2)
        K4 = f(t + h, y + h * K3)
        y = y + (h / 6) * (K1 + 2 * K2 + 2 * K3 + K4)
        t += h
        Y.append(y)

    return Y

def L(sigma, rho, beta):
    """
    Lorenz system
    """
    def f(t, x):
        f1 = sigma * (x[1] - x[0])
        f2 = x[0] * (rho - x[2]) - x[1]
        f3 = x[0] * x[1] - beta * x[2]
        return np.array([f1, f2, f3])
    return f

def test(case): # generate figures in the report.
    x = lambda res: [t[0] for t in res]
    y = lambda res: [t[1] for t in res]
    z = lambda res: [t[2] for t in res]
    h = 0.001
    if case == 1 or case == 2:
        sigma = 10; rho = 28; beta = 8 / 3
    

    if case == 1:
        f = L(sigma, rho, beta)
        y0 = np.array([1.0, 1.0, 1.0])
        y1 = np.array([1 + 1e-4, 1.0, 1.0])
        y2 = np.array([1.0, 1 + 1e-4, 1.0])
        y3 = np.array([1.0, 1.0, 1 + 1e-4])

        fig = plt.figure()

        ax = fig.add_subplot(2, 2, 1, projection='3d')
        T = 10
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(1+1e-3, 1, 1)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(1, 1+1e-3, 1)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(1, 1, 1+1e-3)')
        ax._gridOn
        ax.set_title('T=%d' % T,fontsize=12)
        ax.legend()

        ax = fig.add_subplot(2, 2, 2, projection='3d')
        T = 20
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(1+1e-3, 1, 1)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(1, 1+1e-3, 1)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(1, 1, 1+1e-3)')
        ax._gridOn
        ax.set_title('T=%d' % T,fontsize=12)
        ax.legend()

        ax = fig.add_subplot(2, 2, 3, projection='3d')
        T = 40
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(1+1e-3, 1, 1)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(1, 1+1e-3, 1)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(1, 1, 1+1e-3)')
        ax._gridOn
        ax.set_title('T=%d' % T,fontsize=12)
        ax.legend()

        ax = fig.add_subplot(2, 2, 4, projection='3d')
        T = 60
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(1+1e-3, 1, 1)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(1, 1+1e-3, 1)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(1, 1, 1+1e-3)')
        ax._gridOn
        ax.set_title('T=%d' % T,fontsize=12)
        ax.legend()
        tikzplotlib.save("mytikz1.tex")
        plt.show()
    elif case == 2:
        f = L(sigma, rho, beta)
        y0 = np.array([1.0, 1.0, 1.0])
        y1 = np.array([100.0, 100.0, 100.0])
        y2 = np.array([-100.0, -100.0, -100.0])
        y3 = np.array([-100.0, 100.0, 100.0])

        fig = plt.figure()

        ax = fig.add_subplot(2, 2, 1, projection='3d')
        T = 10
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(100, 100, 100)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(-100, -100, -100)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(-100, 100, 100)')
        ax._gridOn
        ax.set_title('T=%d' % T,fontsize=12)
        ax.legend()

        ax = fig.add_subplot(2, 2, 2, projection='3d')
        T = 20
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(100, 100, 100)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(-100, -100, -100)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(-100, 100, 100)')
        ax._gridOn
        ax.set_title('T=%d' % T,fontsize=12)
        ax.legend()

        ax = fig.add_subplot(2, 2, 3, projection='3d')
        T = 40
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(100, 100, 100)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(-100, -100, -100)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(-100, 100, 100)')
        ax._gridOn
        ax.set_title('T=%d' % T,fontsize=12)
        ax.legend()

        ax = fig.add_subplot(2, 2, 4, projection='3d')
        T = 60
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(100, 100, 100)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(-100, -100, -100)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(-100, 100, 100)')
        ax._gridOn
        ax.set_title('T=%d' % T,fontsize=12)
        ax.legend()
        plt.show()
    elif case == 3:
        sigma = 10; rho = 0.1; beta = 8 / 3
        f = L(sigma, rho, beta)
        y0 = np.array([1.0, 1.0, 1.0])
        y1 = np.array([1.0, 1.0, -1.0])
        y2 = np.array([1.0, -1.0, 1.0])
        y3 = np.array([-1.0, 1.0, 1.0])
        y4 = np.array([1.0, -1.0, -1.0])
        y5 = np.array([-1.0, -1.0, 1.0])
        y6 = np.array([-1.0, 1.0, -1.0])
        y7 = np.array([-1.0, -1.0, -1.0])

        fig = plt.figure()

        ax = fig.add_subplot(2, 2, 1, projection='3d')
        T = 1
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        Y5 = RungeKutta(f, 0, y4, h, T)
        Y6 = RungeKutta(f, 0, y5, h, T)
        Y7 = RungeKutta(f, 0, y6, h, T)
        Y8 = RungeKutta(f, 0, y7, h, T)

        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(1, 1, -1)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(1, -1, 1)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(-1, 1, 1)')
        ax.plot(x(Y5), y(Y5), z(Y5), label = '(1, -1, -1)')
        ax.plot(x(Y6), y(Y6), z(Y6), label = '(-1, -1, 1)')
        ax.plot(x(Y7), y(Y7), z(Y7), label = '(-1, 1, -1)')
        ax.plot(x(Y8), y(Y8), z(Y8), label = '(-1, -1, -1)')

        ax._gridOn
        ax.set_title(r'$\sigma$='+str(sigma)+ r', $\rho$='+str(rho)+r', $\beta=\frac{8}{3}$, '+r'$T$='+str(T),fontsize=12)
        ax.legend()

        ax = fig.add_subplot(2, 2, 2, projection='3d')
        T = 10
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        Y5 = RungeKutta(f, 0, y4, h, T)
        Y6 = RungeKutta(f, 0, y5, h, T)
        Y7 = RungeKutta(f, 0, y6, h, T)
        Y8 = RungeKutta(f, 0, y7, h, T)

        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(1, 1, -1)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(1, -1, 1)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(-1, 1, 1)')
        ax.plot(x(Y5), y(Y5), z(Y5), label = '(1, -1, -1)')
        ax.plot(x(Y6), y(Y6), z(Y6), label = '(-1, -1, 1)')
        ax.plot(x(Y7), y(Y7), z(Y7), label = '(-1, 1, -1)')
        ax.plot(x(Y8), y(Y8), z(Y8), label = '(-1, -1, -1)')
        ax._gridOn
        ax.set_title(r'$\sigma$='+str(sigma)+ r', $\rho$='+str(rho)+r', $\beta=\frac{8}{3}$, '+r'$T$='+str(T),fontsize=12)
        ax.legend()

        sigma = 10; rho = 300; beta = 8 / 3
        f = L(sigma, rho, beta)

        ax = fig.add_subplot(2, 2, 3, projection='3d')
        T = 1
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        Y5 = RungeKutta(f, 0, y4, h, T)
        Y6 = RungeKutta(f, 0, y5, h, T)
        Y7 = RungeKutta(f, 0, y6, h, T)
        Y8 = RungeKutta(f, 0, y7, h, T)

        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(1, 1, -1)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(1, -1, 1)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(-1, 1, 1)')
        ax.plot(x(Y5), y(Y5), z(Y5), label = '(1, -1, -1)')
        ax.plot(x(Y6), y(Y6), z(Y6), label = '(-1, -1, 1)')
        ax.plot(x(Y7), y(Y7), z(Y7), label = '(-1, 1, -1)')
        ax.plot(x(Y8), y(Y8), z(Y8), label = '(-1, -1, -1)')
        ax._gridOn
        ax.set_title(r'$\sigma$='+str(sigma)+ r', $\rho$='+str(rho)+r', $\beta=\frac{8}{3}$, '+r'$T$='+str(T),fontsize=12)
        ax.legend()

        ax = fig.add_subplot(2, 2, 4, projection='3d')
        T = 10
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        Y5 = RungeKutta(f, 0, y4, h, T)
        Y6 = RungeKutta(f, 0, y5, h, T)
        Y7 = RungeKutta(f, 0, y6, h, T)
        Y8 = RungeKutta(f, 0, y7, h, T)

        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(1, 1, -1)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(1, -1, 1)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(-1, 1, 1)')
        ax.plot(x(Y5), y(Y5), z(Y5), label = '(1, -1, -1)')
        ax.plot(x(Y6), y(Y6), z(Y6), label = '(-1, -1, 1)')
        ax.plot(x(Y7), y(Y7), z(Y7), label = '(-1, 1, -1)')
        ax.plot(x(Y8), y(Y8), z(Y8), label = '(-1, -1, -1)')

        ax._gridOn
        ax.set_title(r'$\sigma$='+str(sigma)+ r', $\rho$='+str(rho)+r', $\beta=\frac{8}{3}$, '+r'$T$='+str(T),fontsize=12)
        ax.legend()
        plt.show()
    elif case == 4:


        sigma = 10; rho = 28; beta = 0.8 / 3
        f = L(sigma, rho, beta)
        fig = plt.figure()
        ax = fig.add_subplot(2, 2, 1, projection='3d')
        T = 1
        y0 = np.array([1.0, 1.0, 1.0])
        y1 = np.array([1.0, 1.0, -1.0])
        y2 = np.array([1.0, -1.0, 1.0])
        y3 = np.array([-1.0, 1.0, 1.0])
        y4 = np.array([1.0, -1.0, -1.0])
        y5 = np.array([-1.0, -1.0, 1.0])
        y6 = np.array([-1.0, 1.0, -1.0])
        y7 = np.array([-1.0, -1.0, -1.0])
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        Y5 = RungeKutta(f, 0, y4, h, T)
        Y6 = RungeKutta(f, 0, y5, h, T)
        Y7 = RungeKutta(f, 0, y6, h, T)
        Y8 = RungeKutta(f, 0, y7, h, T)

        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(1, 1, -1)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(1, -1, 1)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(-1, 1, 1)')
        ax.plot(x(Y5), y(Y5), z(Y5), label = '(1, -1, -1)')
        ax.plot(x(Y6), y(Y6), z(Y6), label = '(-1, -1, 1)')
        ax.plot(x(Y7), y(Y7), z(Y7), label = '(-1, 1, -1)')
        ax.plot(x(Y8), y(Y8), z(Y8), label = '(-1, -1, -1)')
        ax._gridOn
        ax.set_title(r'$\sigma$='+str(sigma)+ r', $\rho$='+str(rho)+r', $\beta=\frac{0.8}{3}$, '+r'$T$='+str(T),fontsize=12)
        ax.legend()

        ax = fig.add_subplot(2, 2, 2, projection='3d')
        T = 10
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        Y5 = RungeKutta(f, 0, y4, h, T)
        Y6 = RungeKutta(f, 0, y5, h, T)
        Y7 = RungeKutta(f, 0, y6, h, T)
        Y8 = RungeKutta(f, 0, y7, h, T)

        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(1, 1, -1)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(1, -1, 1)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(-1, 1, 1)')
        ax.plot(x(Y5), y(Y5), z(Y5), label = '(1, -1, -1)')
        ax.plot(x(Y6), y(Y6), z(Y6), label = '(-1, -1, 1)')
        ax.plot(x(Y7), y(Y7), z(Y7), label = '(-1, 1, -1)')
        ax.plot(x(Y8), y(Y8), z(Y8), label = '(-1, -1, -1)')
        
        ax._gridOn
        ax.set_title(r'$\sigma$='+str(sigma)+ r', $\rho$='+str(rho)+r', $\beta=\frac{0.8}{3}$, '+r'$T$='+str(T),fontsize=12)
        ax.legend()

        sigma = 10; rho = 28; beta = 80 / 3
        f = L(sigma, rho, beta)

        ax = fig.add_subplot(2, 2, 3, projection='3d')
        T = 1
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        Y5 = RungeKutta(f, 0, y4, h, T)
        Y6 = RungeKutta(f, 0, y5, h, T)
        Y7 = RungeKutta(f, 0, y6, h, T)
        Y8 = RungeKutta(f, 0, y7, h, T)

        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(1, 1, -1)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(1, -1, 1)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(-1, 1, 1)')
        ax.plot(x(Y5), y(Y5), z(Y5), label = '(1, -1, -1)')
        ax.plot(x(Y6), y(Y6), z(Y6), label = '(-1, -1, 1)')
        ax.plot(x(Y7), y(Y7), z(Y7), label = '(-1, 1, -1)')
        ax.plot(x(Y8), y(Y8), z(Y8), label = '(-1, -1, -1)')
        ax._gridOn
        ax.set_title(r'$\sigma$='+str(sigma)+ r', $\rho$='+str(rho)+r', $\beta=\frac{80}{3}$, '+r'$T$='+str(T),fontsize=12)
        ax.legend()

        ax = fig.add_subplot(2, 2, 4, projection='3d')
        T = 10
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        Y5 = RungeKutta(f, 0, y4, h, T)
        Y6 = RungeKutta(f, 0, y5, h, T)
        Y7 = RungeKutta(f, 0, y6, h, T)
        Y8 = RungeKutta(f, 0, y7, h, T)

        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(1, 1, -1)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(1, -1, 1)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(-1, 1, 1)')
        ax.plot(x(Y5), y(Y5), z(Y5), label = '(1, -1, -1)')
        ax.plot(x(Y6), y(Y6), z(Y6), label = '(-1, -1, 1)')
        ax.plot(x(Y7), y(Y7), z(Y7), label = '(-1, 1, -1)')
        ax.plot(x(Y8), y(Y8), z(Y8), label = '(-1, -1, -1)')
        
        ax._gridOn
        ax.set_title(r'$\sigma$='+str(sigma)+ r', $\rho$='+str(rho)+r', $\beta=\frac{80}{3}$, '+r'$T$='+str(T),fontsize=12)
        ax.legend()

        plt.show()
    elif case == 5:
        sigma = 0.1; rho = 28; beta = 8 / 3
        f = L(sigma, rho, beta)
        y0 = np.array([1.0, 1.0, 1.0])
        y1 = np.array([1.0, 1.0, -1.0])
        y2 = np.array([1.0, -1.0, 1.0])
        y3 = np.array([-1.0, 1.0, 1.0])
        y4 = np.array([1.0, -1.0, -1.0])
        y5 = np.array([-1.0, -1.0, 1.0])
        y6 = np.array([-1.0, 1.0, -1.0])
        y7 = np.array([-1.0, -1.0, -1.0])

        fig = plt.figure()

        ax = fig.add_subplot(3, 2, 1, projection='3d')
        T = 1
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        Y5 = RungeKutta(f, 0, y4, h, T)
        Y6 = RungeKutta(f, 0, y5, h, T)
        Y7 = RungeKutta(f, 0, y6, h, T)
        Y8 = RungeKutta(f, 0, y7, h, T)

        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(1, 1, -1)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(1, -1, 1)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(-1, 1, 1)')
        ax.plot(x(Y5), y(Y5), z(Y5), label = '(1, -1, -1)')
        ax.plot(x(Y6), y(Y6), z(Y6), label = '(-1, -1, 1)')
        ax.plot(x(Y7), y(Y7), z(Y7), label = '(-1, 1, -1)')
        ax.plot(x(Y8), y(Y8), z(Y8), label = '(-1, -1, -1)')

        ax._gridOn
        ax.set_title(r'$\sigma$='+str(sigma)+ r', $\rho$='+str(rho)+r', $\beta=\frac{8}{3}$, '+r'$T$='+str(T),fontsize=12)
        ax.legend()

        ax = fig.add_subplot(3, 2, 2, projection='3d')
        T = 10
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        Y5 = RungeKutta(f, 0, y4, h, T)
        Y6 = RungeKutta(f, 0, y5, h, T)
        Y7 = RungeKutta(f, 0, y6, h, T)
        Y8 = RungeKutta(f, 0, y7, h, T)

        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(1, 1, -1)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(1, -1, 1)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(-1, 1, 1)')
        ax.plot(x(Y5), y(Y5), z(Y5), label = '(1, -1, -1)')
        ax.plot(x(Y6), y(Y6), z(Y6), label = '(-1, -1, 1)')
        ax.plot(x(Y7), y(Y7), z(Y7), label = '(-1, 1, -1)')
        ax.plot(x(Y8), y(Y8), z(Y8), label = '(-1, -1, -1)')
        ax._gridOn
        ax.set_title(r'$\sigma$='+str(sigma)+ r', $\rho$='+str(rho)+r', $\beta=\frac{8}{3}$, '+r'$T$='+str(T),fontsize=12)
        ax.legend()

        sigma = 1; rho = 28; beta = 8 / 3
        f = L(sigma, rho, beta)

        ax = fig.add_subplot(3, 2, 3, projection='3d')
        T = 1
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        Y5 = RungeKutta(f, 0, y4, h, T)
        Y6 = RungeKutta(f, 0, y5, h, T)
        Y7 = RungeKutta(f, 0, y6, h, T)
        Y8 = RungeKutta(f, 0, y7, h, T)

        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(1, 1, -1)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(1, -1, 1)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(-1, 1, 1)')
        ax.plot(x(Y5), y(Y5), z(Y5), label = '(1, -1, -1)')
        ax.plot(x(Y6), y(Y6), z(Y6), label = '(-1, -1, 1)')
        ax.plot(x(Y7), y(Y7), z(Y7), label = '(-1, 1, -1)')
        ax.plot(x(Y8), y(Y8), z(Y8), label = '(-1, -1, -1)')

        ax._gridOn
        ax.set_title(r'$\sigma$='+str(sigma)+ r', $\rho$='+str(rho)+r', $\beta=\frac{8}{3}$, '+r'$T$='+str(T),fontsize=12)
        ax.legend()

        ax = fig.add_subplot(3, 2, 4, projection='3d')
        T = 10
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        Y5 = RungeKutta(f, 0, y4, h, T)
        Y6 = RungeKutta(f, 0, y5, h, T)
        Y7 = RungeKutta(f, 0, y6, h, T)
        Y8 = RungeKutta(f, 0, y7, h, T)

        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(1, 1, -1)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(1, -1, 1)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(-1, 1, 1)')
        ax.plot(x(Y5), y(Y5), z(Y5), label = '(1, -1, -1)')
        ax.plot(x(Y6), y(Y6), z(Y6), label = '(-1, -1, 1)')
        ax.plot(x(Y7), y(Y7), z(Y7), label = '(-1, 1, -1)')
        ax.plot(x(Y8), y(Y8), z(Y8), label = '(-1, -1, -1)')
        ax._gridOn
        ax.set_title(r'$\sigma$='+str(sigma)+ r', $\rho$='+str(rho)+r', $\beta=\frac{8}{3}$, '+r'$T$='+str(T),fontsize=12)
        ax.legend()

        sigma = 100; rho = 28; beta = 8 / 3
        f = L(sigma, rho, beta)

        ax = fig.add_subplot(3, 2, 5, projection='3d')
        T = 1
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        Y5 = RungeKutta(f, 0, y4, h, T)
        Y6 = RungeKutta(f, 0, y5, h, T)
        Y7 = RungeKutta(f, 0, y6, h, T)
        Y8 = RungeKutta(f, 0, y7, h, T)

        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(1, 1, -1)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(1, -1, 1)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(-1, 1, 1)')
        ax.plot(x(Y5), y(Y5), z(Y5), label = '(1, -1, -1)')
        ax.plot(x(Y6), y(Y6), z(Y6), label = '(-1, -1, 1)')
        ax.plot(x(Y7), y(Y7), z(Y7), label = '(-1, 1, -1)')
        ax.plot(x(Y8), y(Y8), z(Y8), label = '(-1, -1, -1)')
        ax._gridOn
        ax.set_title(r'$\sigma$='+str(sigma)+ r', $\rho$='+str(rho)+r', $\beta=\frac{8}{3}$, '+r'$T$='+str(T),fontsize=12)
        ax.legend()

        ax = fig.add_subplot(3, 2, 6, projection='3d')
        T = 10
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        Y5 = RungeKutta(f, 0, y4, h, T)
        Y6 = RungeKutta(f, 0, y5, h, T)
        Y7 = RungeKutta(f, 0, y6, h, T)
        Y8 = RungeKutta(f, 0, y7, h, T)

        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(1, 1, -1)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(1, -1, 1)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(-1, 1, 1)')
        ax.plot(x(Y5), y(Y5), z(Y5), label = '(1, -1, -1)')
        ax.plot(x(Y6), y(Y6), z(Y6), label = '(-1, -1, 1)')
        ax.plot(x(Y7), y(Y7), z(Y7), label = '(-1, 1, -1)')
        ax.plot(x(Y8), y(Y8), z(Y8), label = '(-1, -1, -1)')

        ax._gridOn
        ax.set_title(r'$\sigma$='+str(sigma)+ r', $\rho$='+str(rho)+r', $\beta=\frac{8}{3}$, '+r'$T$='+str(T),fontsize=12)
        ax.legend()
        plt.show()
    elif case == 6:
        sigma = 10; rho = 28; beta = 80 / 3
        f = L(sigma, rho, beta)
        y0 = np.array([1.0, 1.0, 1.0])
        y1 = np.array([1.0, 1.0, -1.0])
        y2 = np.array([1.0, -1.0, 1.0])
        y3 = np.array([-1.0, 1.0, 1.0])
        y4 = np.array([1.0, -1.0, -1.0])
        y5 = np.array([-1.0, -1.0, 1.0])
        y6 = np.array([-1.0, 1.0, -1.0])
        y7 = np.array([-1.0, -1.0, -1.0])

        fig = plt.figure()

        ax = fig.add_subplot(3, 2, 1, projection='3d')
        T = 1
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        Y5 = RungeKutta(f, 0, y4, h, T)
        Y6 = RungeKutta(f, 0, y5, h, T)
        Y7 = RungeKutta(f, 0, y6, h, T)
        Y8 = RungeKutta(f, 0, y7, h, T)

        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(1, 1, -1)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(1, -1, 1)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(-1, 1, 1)')
        ax.plot(x(Y5), y(Y5), z(Y5), label = '(1, -1, -1)')
        ax.plot(x(Y6), y(Y6), z(Y6), label = '(-1, -1, 1)')
        ax.plot(x(Y7), y(Y7), z(Y7), label = '(-1, 1, -1)')
        ax.plot(x(Y8), y(Y8), z(Y8), label = '(-1, -1, -1)')

        ax._gridOn
        ax.set_title(r'$\sigma$='+str(sigma)+ r', $\rho$='+str(rho)+r', $\beta=\frac{80}{3}$, '+r'$T$='+str(T),fontsize=12)
        ax.legend()

        ax = fig.add_subplot(3, 2, 2, projection='3d')
        T = 10
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        Y5 = RungeKutta(f, 0, y4, h, T)
        Y6 = RungeKutta(f, 0, y5, h, T)
        Y7 = RungeKutta(f, 0, y6, h, T)
        Y8 = RungeKutta(f, 0, y7, h, T)

        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(1, 1, -1)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(1, -1, 1)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(-1, 1, 1)')
        ax.plot(x(Y5), y(Y5), z(Y5), label = '(1, -1, -1)')
        ax.plot(x(Y6), y(Y6), z(Y6), label = '(-1, -1, 1)')
        ax.plot(x(Y7), y(Y7), z(Y7), label = '(-1, 1, -1)')
        ax.plot(x(Y8), y(Y8), z(Y8), label = '(-1, -1, -1)')
        ax._gridOn
        ax.set_title(r'$\sigma$='+str(sigma)+ r', $\rho$='+str(rho)+r', $\beta=\frac{80}{3}$, '+r'$T$='+str(T),fontsize=12)
        ax.legend()

        sigma = 10; rho = 28; beta = 3.2 / 3
        f = L(sigma, rho, beta)

        ax = fig.add_subplot(3, 2, 3, projection='3d')
        T = 1
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        Y5 = RungeKutta(f, 0, y4, h, T)
        Y6 = RungeKutta(f, 0, y5, h, T)
        Y7 = RungeKutta(f, 0, y6, h, T)
        Y8 = RungeKutta(f, 0, y7, h, T)

        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(1, 1, -1)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(1, -1, 1)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(-1, 1, 1)')
        ax.plot(x(Y5), y(Y5), z(Y5), label = '(1, -1, -1)')
        ax.plot(x(Y6), y(Y6), z(Y6), label = '(-1, -1, 1)')
        ax.plot(x(Y7), y(Y7), z(Y7), label = '(-1, 1, -1)')
        ax.plot(x(Y8), y(Y8), z(Y8), label = '(-1, -1, -1)')

        ax._gridOn
        ax.set_title(r'$\sigma$='+str(sigma)+ r', $\rho$='+str(rho)+r', $\beta=\frac{3.2}{3}$, '+r'$T$='+str(T),fontsize=12)
        ax.legend()

        ax = fig.add_subplot(3, 2, 4, projection='3d')
        T = 10
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        Y5 = RungeKutta(f, 0, y4, h, T)
        Y6 = RungeKutta(f, 0, y5, h, T)
        Y7 = RungeKutta(f, 0, y6, h, T)
        Y8 = RungeKutta(f, 0, y7, h, T)

        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(1, 1, -1)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(1, -1, 1)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(-1, 1, 1)')
        ax.plot(x(Y5), y(Y5), z(Y5), label = '(1, -1, -1)')
        ax.plot(x(Y6), y(Y6), z(Y6), label = '(-1, -1, 1)')
        ax.plot(x(Y7), y(Y7), z(Y7), label = '(-1, 1, -1)')
        ax.plot(x(Y8), y(Y8), z(Y8), label = '(-1, -1, -1)')
        ax._gridOn
        ax.set_title(r'$\sigma$='+str(sigma)+ r', $\rho$='+str(rho)+r', $\beta=\frac{3.2}{3}$, '+r'$T$='+str(T),fontsize=12)
        ax.legend()

        sigma = 10; rho = 28; beta = 0.8 / 3
        f = L(sigma, rho, beta)

        ax = fig.add_subplot(3, 2, 5, projection='3d')
        T = 4
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        Y5 = RungeKutta(f, 0, y4, h, T)
        Y6 = RungeKutta(f, 0, y5, h, T)
        Y7 = RungeKutta(f, 0, y6, h, T)
        Y8 = RungeKutta(f, 0, y7, h, T)

        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(1, 1, -1)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(1, -1, 1)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(-1, 1, 1)')
        ax.plot(x(Y5), y(Y5), z(Y5), label = '(1, -1, -1)')
        ax.plot(x(Y6), y(Y6), z(Y6), label = '(-1, -1, 1)')
        ax.plot(x(Y7), y(Y7), z(Y7), label = '(-1, 1, -1)')
        ax.plot(x(Y8), y(Y8), z(Y8), label = '(-1, -1, -1)')
        ax._gridOn
        ax.set_title(r'$\sigma$='+str(sigma)+ r', $\rho$='+str(rho)+r', $\beta=\frac{0.8}{3}$, '+r'$T$='+str(T),fontsize=12)
        ax.legend()

        ax = fig.add_subplot(3, 2, 6, projection='3d')
        T = 80
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        Y5 = RungeKutta(f, 0, y4, h, T)
        Y6 = RungeKutta(f, 0, y5, h, T)
        Y7 = RungeKutta(f, 0, y6, h, T)
        Y8 = RungeKutta(f, 0, y7, h, T)

        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(1, 1, -1)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(1, -1, 1)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(-1, 1, 1)')
        ax.plot(x(Y5), y(Y5), z(Y5), label = '(1, -1, -1)')
        ax.plot(x(Y6), y(Y6), z(Y6), label = '(-1, -1, 1)')
        ax.plot(x(Y7), y(Y7), z(Y7), label = '(-1, 1, -1)')
        ax.plot(x(Y8), y(Y8), z(Y8), label = '(-1, -1, -1)')

        ax._gridOn
        ax.set_title(r'$\sigma$='+str(sigma)+ r', $\rho$='+str(rho)+r', $\beta=\frac{0.8}{3}$, '+r'$T$='+str(T),fontsize=12)
        ax.legend()
        plt.show()

    elif case == 7:


        sigma = -1; rho = 28; beta = 8 / 3
        f = L(sigma, rho, beta)
        fig = plt.figure()
        ax = fig.add_subplot(3, 2, 1, projection='3d')
        T = 1
        y0 = np.array([1.0, 1.0, 1.0])
        y1 = np.array([1.0, 1.0, -1.0])
        y2 = np.array([1.0, -1.0, 1.0])
        y3 = np.array([-1.0, 1.0, 1.0])
        y4 = np.array([1.0, -1.0, -1.0])
        y5 = np.array([-1.0, -1.0, 1.0])
        y6 = np.array([-1.0, 1.0, -1.0])
        y7 = np.array([-1.0, -1.0, -1.0])
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        Y5 = RungeKutta(f, 0, y4, h, T)
        Y6 = RungeKutta(f, 0, y5, h, T)
        Y7 = RungeKutta(f, 0, y6, h, T)
        Y8 = RungeKutta(f, 0, y7, h, T)

        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(1, 1, -1)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(1, -1, 1)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(-1, 1, 1)')
        ax.plot(x(Y5), y(Y5), z(Y5), label = '(1, -1, -1)')
        ax.plot(x(Y6), y(Y6), z(Y6), label = '(-1, -1, 1)')
        ax.plot(x(Y7), y(Y7), z(Y7), label = '(-1, 1, -1)')
        ax.plot(x(Y8), y(Y8), z(Y8), label = '(-1, -1, -1)')
        ax._gridOn
        ax.set_title(r'$\sigma$='+str(sigma)+ r', $\rho$='+str(rho)+r', $\beta=\frac{8}{3}$, '+r'$T$='+str(T),fontsize=12)
        ax.legend()

        ax = fig.add_subplot(3, 2, 2, projection='3d')
        T = 10
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        Y5 = RungeKutta(f, 0, y4, h, T)
        Y6 = RungeKutta(f, 0, y5, h, T)
        Y7 = RungeKutta(f, 0, y6, h, T)
        Y8 = RungeKutta(f, 0, y7, h, T)

        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(1, 1, -1)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(1, -1, 1)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(-1, 1, 1)')
        ax.plot(x(Y5), y(Y5), z(Y5), label = '(1, -1, -1)')
        ax.plot(x(Y6), y(Y6), z(Y6), label = '(-1, -1, 1)')
        ax.plot(x(Y7), y(Y7), z(Y7), label = '(-1, 1, -1)')
        ax.plot(x(Y8), y(Y8), z(Y8), label = '(-1, -1, -1)')
        
        ax._gridOn
        ax.set_title(r'$\sigma$='+str(sigma)+ r', $\rho$='+str(rho)+r', $\beta=\frac{8}{3}$, '+r'$T$='+str(T),fontsize=12)
        ax.legend()

        sigma = 10; rho = -1; beta = 8 / 3
        f = L(sigma, rho, beta)

        ax = fig.add_subplot(3, 2, 3, projection='3d')
        T = 1
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        Y5 = RungeKutta(f, 0, y4, h, T)
        Y6 = RungeKutta(f, 0, y5, h, T)
        Y7 = RungeKutta(f, 0, y6, h, T)
        Y8 = RungeKutta(f, 0, y7, h, T)

        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(1, 1, -1)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(1, -1, 1)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(-1, 1, 1)')
        ax.plot(x(Y5), y(Y5), z(Y5), label = '(1, -1, -1)')
        ax.plot(x(Y6), y(Y6), z(Y6), label = '(-1, -1, 1)')
        ax.plot(x(Y7), y(Y7), z(Y7), label = '(-1, 1, -1)')
        ax.plot(x(Y8), y(Y8), z(Y8), label = '(-1, -1, -1)')
        ax._gridOn
        ax.set_title(r'$\sigma$='+str(sigma)+ r', $\rho$='+str(rho)+r', $\beta=\frac{8}{3}$, '+r'$T$='+str(T),fontsize=12)
        ax.legend()

        ax = fig.add_subplot(3, 2, 4, projection='3d')
        T = 10
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        Y5 = RungeKutta(f, 0, y4, h, T)
        Y6 = RungeKutta(f, 0, y5, h, T)
        Y7 = RungeKutta(f, 0, y6, h, T)
        Y8 = RungeKutta(f, 0, y7, h, T)

        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(1, 1, -1)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(1, -1, 1)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(-1, 1, 1)')
        ax.plot(x(Y5), y(Y5), z(Y5), label = '(1, -1, -1)')
        ax.plot(x(Y6), y(Y6), z(Y6), label = '(-1, -1, 1)')
        ax.plot(x(Y7), y(Y7), z(Y7), label = '(-1, 1, -1)')
        ax.plot(x(Y8), y(Y8), z(Y8), label = '(-1, -1, -1)')
        
        ax._gridOn
        ax.set_title(r'$\sigma$='+str(sigma)+ r', $\rho$='+str(rho)+r', $\beta=\frac{8}{3}$, '+r'$T$='+str(T),fontsize=12)
        ax.legend()

        sigma = 10; rho = 28; beta = -8 / 3
        f = L(sigma, rho, beta)

        ax = fig.add_subplot(3, 2, 5, projection='3d')
        T = 1
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        Y5 = RungeKutta(f, 0, y4, h, T)
        Y6 = RungeKutta(f, 0, y5, h, T)
        Y7 = RungeKutta(f, 0, y6, h, T)
        Y8 = RungeKutta(f, 0, y7, h, T)

        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(1, 1, -1)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(1, -1, 1)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(-1, 1, 1)')
        ax.plot(x(Y5), y(Y5), z(Y5), label = '(1, -1, -1)')
        ax.plot(x(Y6), y(Y6), z(Y6), label = '(-1, -1, 1)')
        ax.plot(x(Y7), y(Y7), z(Y7), label = '(-1, 1, -1)')
        ax.plot(x(Y8), y(Y8), z(Y8), label = '(-1, -1, -1)')
        ax._gridOn
        ax.set_title(r'$\sigma$='+str(sigma)+ r', $\rho$='+str(rho)+r', $\beta=-\frac{8}{3}$, '+r'$T$='+str(T),fontsize=12)
        ax.legend()

        ax = fig.add_subplot(3, 2, 6, projection='3d')
        T = 10
        Y1 = RungeKutta(f, 0, y0, h, T)
        Y2 = RungeKutta(f, 0, y1, h, T)
        Y3 = RungeKutta(f, 0, y2, h, T)
        Y4 = RungeKutta(f, 0, y3, h, T)
        Y5 = RungeKutta(f, 0, y4, h, T)
        Y6 = RungeKutta(f, 0, y5, h, T)
        Y7 = RungeKutta(f, 0, y6, h, T)
        Y8 = RungeKutta(f, 0, y7, h, T)

        ax.plot(x(Y1), y(Y1), z(Y1), label = '(1, 1, 1)')
        ax.plot(x(Y2), y(Y2), z(Y2), label = '(1, 1, -1)')
        ax.plot(x(Y3), y(Y3), z(Y3), label = '(1, -1, 1)')
        ax.plot(x(Y4), y(Y4), z(Y4), label = '(-1, 1, 1)')
        ax.plot(x(Y5), y(Y5), z(Y5), label = '(1, -1, -1)')
        ax.plot(x(Y6), y(Y6), z(Y6), label = '(-1, -1, 1)')
        ax.plot(x(Y7), y(Y7), z(Y7), label = '(-1, 1, -1)')
        ax.plot(x(Y8), y(Y8), z(Y8), label = '(-1, -1, -1)')
        
        ax._gridOn
        ax.set_title(r'$\sigma$='+str(sigma)+ r', $\rho$='+str(rho)+r', $\beta=-\frac{8}{3}$, '+r'$T$='+str(T),fontsize=12)
        ax.legend()

        plt.show()



if __name__ == "__main__":
    for i in range(1,8):
        test(i)