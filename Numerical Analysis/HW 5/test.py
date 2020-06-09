from ode_solver import *
from schemes_adp import *
import numpy as np

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

def compare():
    sigma = 10; rho = 28; beta = 8 / 3
    f = L(sigma, rho, beta)
    y0 = np.array([1.0, 1.0, 1.0])
    T = 1
    Err = []
    h = 0.0001
    h0 = 0.0001 # 0.0001 T=1,10; 0.001 T=100
    Y0 = ode45(f, 0, y0, h0, T); Y0 = Y0[0]
    Y01 = ode45(f, 0, y0, h0, T, 'Fehlberg'); Y01 = Y01[0]; print(np.linalg.norm(Y01[-1]-Y0[-1]))
    Y02 = ode45(f, 0, y0, h0, T, 'Cash_Karp'); Y02 = Y02[0]; print(np.linalg.norm(Y02[-1]-Y0[-1]))
    
    # runge-kutta
    
    Y = ode_solver(f, 0, y0, h, T,'Kutta33'); Err.append(np.linalg.norm(Y[-1]-Y0[-1]))
    Y = ode_solver(f, 0, y0, h, T,'Heun33'); Err.append(np.linalg.norm(Y[-1]-Y0[-1]))
    Y = ode_solver(f, 0, y0, h, T,'RungeKutta'); Err.append(np.linalg.norm(Y[-1]-Y0[-1]))
    Y = ode_solver(f, 0, y0, h, T,'Kutta44'); Err.append(np.linalg.norm(Y[-1]-Y0[-1]))
    Y = ode_solver(f, 0, y0, h, T,'Gill44'); Err.append(np.linalg.norm(Y[-1]-Y0[-1]))
    Y = ode_solver(f, 0, y0, h, T,'Bogacki_Shampine_3'); Err.append(np.linalg.norm(Y[-1]-Y0[-1]))
    Y = ode_solver(f, 0, y0, h, T,'Bogacki_Shampine_2'); Err.append(np.linalg.norm(Y[-1]-Y0[-1]))
    Y = ode_solver(f, 0, y0, h, T,'Fehlberg_5'); Err.append(np.linalg.norm(Y[-1]-Y0[-1]))
    Y = ode_solver(f, 0, y0, h, T,'Fehlberg_4'); Err.append(np.linalg.norm(Y[-1]-Y0[-1]))
    Y = ode_solver(f, 0, y0, h, T,'Cash_Karp_5'); Err.append(np.linalg.norm(Y[-1]-Y0[-1]))
    Y = ode_solver(f, 0, y0, h, T, 'Cash_Karp_4'); Err.append(np.linalg.norm(Y[-1]-Y0[-1]))
    Y = ode_solver(f, 0, y0, h, T,'Dormand_Prince_5'); Err.append(np.linalg.norm(Y[-1]-Y0[-1]))
    Y = ode_solver(f, 0, y0, h, T, 'Dormand_Prince_4'); Err.append(np.linalg.norm(Y[-1]-Y0[-1]))
    # forward scheme
    Y = Euler(f, 0, y0, h, T); Err.append(np.linalg.norm(Y[-1]-Y0[-1]))
    Y = Adams2(f, 0, y0, h, T); Err.append(np.linalg.norm(Y[-1]-Y0[-1]))
    
    # predict_correct
    Y = Euler2(f, 0, y0, h, T); Err.append(np.linalg.norm(Y[-1]-Y0[-1]))
    Y = Milne3(f, 0, y0, h, T); Err.append(np.linalg.norm(Y[-1]-Y0[-1]))
    Y = Adams4(f, 0, y0, h, T); Err.append(np.linalg.norm(Y[-1]-Y0[-1]))

    return Err

x = compare()