# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 14:19:17 2020

@author: DELL
"""

from schemes import *
def ode_solver(f, t, y0, h, T, method='RungeKutta'):
    method = method + '_step'
    ode_step = eval(method)
    Y = [y0]; ts = [t]
    while t < T:
        Y = ode_step(f, t, Y, h)
        t += h
        ts.append(t)
    Y.pop();
    t = ts[-2]
    h = T - t
    Y = ode_step(f, t, Y, h)
    return Y
