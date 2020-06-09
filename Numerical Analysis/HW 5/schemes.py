from math import sqrt


def RungeKutta_step(f, t, Y, h):
    y = Y[-1]
    K1 = f(t, y)
    K2 = f(t + h / 2, y + (h / 2) * K1)
    K3 = f(t + h / 2, y + (h / 2) * K2)
    K4 = f(t + h, y + h * K3)
    y = y + (h / 6) * (K1 + 2 * K2 + 2 * K3 + K4)
    Y.append(y)
    return Y

def Kutta33_step(f, t, Y, h):
    y = Y[-1]
    K1 = f(t, y)
    K2 = f(t + h / 2, y + (h / 2) * K1)
    K3 = f(t + h, y - h * K1 + 2 * h * K2)
    y = y + (h / 6) * (K1 + 4 * K2 +  K3)
    Y.append(y)
    return Y

def Heun33_step(f, t, Y, h):
    y = Y[-1]
    K1 = f(t, y)
    K2 = f(t + h / 3, y + (h / 3) * K1)
    K3 = f(t + 2 * h / 3, y + (2 * h / 3) * K2)
    y = y + (h / 4) * (K1 + 3 * K3)
    Y.append(y)
    return Y

def Kutta44_step(f, t, Y, h):
    y = Y[-1]
    K1 = f(t, y)
    K2 = f(t + h / 3, y + (h / 3) * K1)
    K3 = f(t + 2 * h / 3, y - (h / 3) * K1 + h * K2)
    K4 = f(t + h, y + h * (K1 - K2 + K3))
    y = y + (h / 8) * (K1 + 3 * K2 + 3 * K3 + K4)
    Y.append(y)
    return Y

def Gill44_step(f, t, Y, h):
    y = Y[-1]
    K1 = f(t, y)
    K2 = f(t + h / 2, y + (h / 2) * K1)
    K3 = f(t + h / 2, y + (sqrt(2)-1)/2 * h * K1 + (1-sqrt(2)/2) * h * K2)
    K4 = f(t + h, y - sqrt(2)/2 * h * K2 + (1+sqrt(2)/2) * h * K3)
    y = y + (h / 6) * (K1 + (2-sqrt(2)) * K2 + (2+sqrt(2)) * K3 + K4)
    Y.append(y)
    return Y


def Euler_step(f, t, Y, h):
    y = Y[-1]
    y = y + h * f(t, y)
    t += h
    Y.append(y)
    return Y

def Euler2_step(f, t, Y, h):
    y = Y[-1]
    y0 = y + h * f(t, y)
    y = y + h / 2 * (f(t, y) + f(t + h, y0))
    Y.append(y)
    return Y


def Bogacki_Shampine_3_step(f, t, Y, h):
    y = Y[-1]
    K1 = f(t, y)
    K2 = f(t + h / 2, y + h / 2 * K1)
    K3 = f(t + 3 * h / 4, y + (h * 3 / 4) * K2)
    K4 = f(t + h, y + (h * 2 / 9) * K1 + (h * 1 / 3) * K2 + (h * 4 / 9) * K3)
    y = y + h * (K1 * 2 / 9 + 1 / 3 * K2 + 4 / 9 * K3)
    Y.append(y)
    return Y

def Bogacki_Shampine_2_step(f, t, Y, h):
    y = Y[-1]
    K1 = f(t, y)
    K2 = f(t + h / 2, y + h / 2 * K1)
    K3 = f(t + 3 * h / 4, y + (h * 3 / 4) * K2)
    K4 = f(t + h, y + (h * 2 / 9) * K1 + (h * 1 / 3) * K2 + (h * 4 / 9) * K3)
    y = y + h * (K1 * 7 / 24 + 1 / 4 * K2 + 1 / 3 * K3 + 1 / 8 * K4)
    Y.append(y)
    return Y

def Fehlberg_5_step(f, t, Y, h):
    y = Y[-1]
    K1 = f(t, y)
    K2 = f(t + h / 4, y + h / 4 * K1)
    K3 = f(t + 3 * h / 8, y + (h * 3 / 32) * K1 + (h * 9 / 32) * K2)
    K4 = f(t + 12 * h / 13, y + (h * 1932 / 2197) * K1 - (h * 7200 / 2197) * K2 + (h * 7296 / 2197) * K3)
    K5 = f(t + h, y + (h * 439 / 216) * K1 - (h * 8) * K2 + (h * 3680 / 513) * K3 - (h * 845 / 4104) * K4)
    K6 = f(t + h / 2, y - (h * 8 / 27) * K1 + (h * 2) * K2 - (h * 3544 / 2565) * K3 + (h * 1859 / 4104) * K4 - (h * 11 / 40) * K5)
    y = y + h * (K1 * 16 / 135 + 6656 / 12825 * K3 + 28561 / 56430 * K4 - 9 / 50 * K5 + 2 / 55 * K6)
    Y.append(y)
    return Y

def Fehlberg_4_step(f, t, Y, h):
    y = Y[-1]
    K1 = f(t, y)
    K2 = f(t + h / 4, y + h / 4 * K1)
    K3 = f(t + 3 * h / 8, y + (h * 3 / 32) * K1 + (h * 9 / 32) * K2)
    K4 = f(t + 12 * h / 13, y + (h * 1932 / 2197) * K1 - (h * 7200 / 2197) * K2 + (h * 7296 / 2197) * K3)
    K5 = f(t + h, y + (h * 439 / 216) * K1 - (h * 8) * K2 + (h * 3680 / 513) * K3 - (h * 845 / 4104) * K4)
    K6 = f(t + h / 2, y - (h * 8 / 27) * K1 + (h * 2) * K2 - (h * 3544 / 2565) * K3 + (h * 1859 / 4104) * K4 - (h * 11 / 40) * K5)
    y = y + h * (K1 * 25 / 216 + 1408 / 2565 * K3 + 2197 / 4104 * K4 - 1 / 5 * K5)
    Y.append(y)
    return Y

def Cash_Karp_5_step(f, t, Y, h):
    y = Y[-1]
    K1 = f(t, y)
    K2 = f(t + h / 5, y + h / 5 * K1)
    K3 = f(t + 3 * h / 10, y + (h * 3 / 40) * K1 + (h * 9 / 40) * K2)
    K4 = f(t + 3 * h / 5, y + (h * 3 / 10) * K1 - (h * 9 / 10) * K2 + (h * 6 / 5) * K3)
    K5 = f(t + h, y - (h * 11 / 54) * K1 + (h * 5 / 2) * K2 - (h * 70 / 27) * K3 + (h * 35 / 27) * K4)
    K6 = f(t + 7 * h / 8, y + (h * 1631 / 55296) * K1 + (h * 175 / 512) * K2 + (h * 575 / 13824) * K3 + (h * 44275 / 110592) * K4 + (h * 253 / 4096) * K5)
    y = y + h * (K1 * 37 / 378 + 250 / 621 * K3 + 125 / 594 * K4 +  512 / 1771 * K6)
    Y.append(y)
    return Y

def Cash_Karp_4_step(f, t, Y, h):
    y = Y[-1]
    K1 = f(t, y)
    K2 = f(t + h / 5, y + h / 5 * K1)
    K3 = f(t + 3 * h / 10, y + (h * 3 / 40) * K1 + (h * 9 / 40) * K2)
    K4 = f(t + 3 * h / 5, y + (h * 3 / 10) * K1 - (h * 9 / 10) * K2 + (h * 6 / 5) * K3)
    K5 = f(t + h, y - (h * 11 / 54) * K1 + (h * 5 / 2) * K2 - (h * 70 / 27) * K3 + (h * 35 / 27) * K4)
    K6 = f(t + 7 * h / 8, y + (h * 1631 / 55296) * K1 + (h * 175 / 512) * K2 + (h * 575 / 13824) * K3 + (h * 44275 / 110592) * K4 + (h * 253 / 4096) * K5)
    y = y + h * (K1 * 2825 / 27648 + 18575 / 48384 * K3 + 13525 / 55296 * K4 + 277 / 14336 * K5 + 1 / 4 * K6)
    Y.append(y)
    return Y
        
def Dormand_Prince_5_step(f, t, Y, h):
    y = Y[-1]
    K1 = f(t, y)
    K2 = f(t + h / 5, y + h / 5 * K1)
    K3 = f(t + 3 * h / 10, y + (h * 3 / 40) * K1 + (h * 9 / 40) * K2)
    K4 = f(t + 4 * h / 5, y + (h * 44 / 45) * K1 - (h * 56 / 15) * K2 + (h * 32 / 9) * K3)
    K5 = f(t + 8 / 9 * h, y + (h * 19372 / 6561) * K1 - (h * 25360 / 2187) * K2 + (h * 64448 / 6561) * K3 - (h * 212 / 729) * K4)
    K6 = f(t + h, y + (h * 9017 / 3168) * K1 - (h * 355 / 33) * K2 + (h * 46732 / 5247) * K3 + (h * 49 / 176) * K4 - (h * 5103 / 18656) * K5)
    K7 = f(t + h, y + (h * 35 / 384) * K1 + (h * 500 / 1113) * K3 + (h * 125 / 192) * K4 - (h * 2187 / 6784) * K5 + (h * 11 / 84) * K6)
    y = y + h * (K1 * 35 / 384 + 500 / 1113 * K3 + 125 / 192 * K4  -  2187 / 6784 * K5 + 11 / 84 * K6)
    Y.append(y)
    return Y

def Dormand_Prince_4_step(f, t, Y, h):
    y = Y[-1]
    K1 = f(t, y)
    K2 = f(t + h / 5, y + h / 5 * K1)
    K3 = f(t + 3 * h / 10, y + (h * 3 / 40) * K1 + (h * 9 / 40) * K2)
    K4 = f(t + 4 * h / 5, y + (h * 44 / 45) * K1 - (h * 56 / 15) * K2 + (h * 32 / 9) * K3)
    K5 = f(t + 8 / 9 * h, y + (h * 19372 / 6561) * K1 - (h * 25360 / 2187) * K2 + (h * 64448 / 6561) * K3 - (h * 212 / 729) * K4)
    K6 = f(t + h, y + (h * 9017 / 3168) * K1 - (h * 355 / 33) * K2 + (h * 46732 / 5247) * K3 + (h * 49 / 176) * K4 - (h * 5103 / 18656) * K5)
    K7 = f(t + h, y + (h * 35 / 384) * K1 + (h * 500 / 1113) * K3 + (h * 125 / 192) * K4 - (h * 2187 / 6784) * K5 + (h * 11 / 84) * K6)
    y = y + h * (K1 * 5179 / 57600 + 7571 / 16695 * K3 + 393 / 640 * K4  -  92097 / 339200 * K5 + 187 / 2100 * K6 + 1 / 40 * K7)
    Y.append(y)
    return Y


def Multi33_step(f, t, Y, h):
    y = Y[-1]
    if len(Y) <= 1:
        y0 = y + h * f(t, y)
        y = y + h / 2 * (f(t, y) + f(t + h, y0))
    else:
        y = Y[-2] + h / 3 * (7*f(t, y) - 2*f(t-h, Y[-2]) + f(t-2*h, Y[-3]))
    t += h
    Y.append(y)
    return Y

def Adams2_step(f, t, Y, h):
    y = Y[-1]
    if len(Y) <= 1:
        y0 = y + h * f(t, y)
        y = y + h / 2 * (f(t, y) + f(t + h, y0))
    else:
        y = y + h / 2 * (3*f(t, y) - f(t-h, Y[-2]))
    Y.append(y)
    return Y

def Milne3_step(f, t, Y, h):
    y = Y[-1]
    if len(Y) <= 3:
        y0 = y + h * f(t, y)
        y = y + h / 2 * (f(t, y) + f(t + h, y0))
    else:
        y0 = Y[-4] + h / 3 * (8 * f(t, y) - 4 * f(t-h, Y[-2]) + 8 * f(t-2*h, Y[-3]))
        y = Y[-2] + h / 3 * (f(t+h, y0) + 4*f(t,y) + f(t-h,Y[-2]))
    Y.append(y)
    return Y

def Adams4_step(f, t, Y, h):
    y = Y[-1]
    if len(Y) <= 2:
        y0 = y + h * f(t, y)
        y = y + h / 2 * (f(t, y) + f(t + h, y0))
    elif len(Y) == 3:
        y0 = Y[-4] + h / 3 * (8 * f(t, y) - 4 * f(t-h, Y[-2]) + 8 * f(t-2*h, Y[-3]))
        y = Y[-2] + h / 3 * (f(t+h, y0) + 4*f(t,y) + f(t-h,Y[-2]))
    else:
        y0 = y + h / 24 * (55 * f(t, y) - 59 * f(t-h, Y[-2]) + 37 * f(t-2*h, Y[-3]) - 9 * f(t-3*h, Y[-4]))
        y = y + h / 24 * (9 *f(t+h, y0) + 19 * f(t, y) - 5 * f(t-h, Y[-2]) + f(t-2*h, Y[-3]))
    Y.append(y)
    return Y
              
        


