from itertools import cycle
from collections import defaultdict

import numpy as np
from numpy import linalg as LA
from log_reg.commons import Oracles, compute_error
from log_reg.algorithms import GD, GDstr
from log_reg.operators import norm1, norm2sq
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import sklearn.linear_model


def main():
    training, testing = np.load('dataset/training.npz'), np.load(
        'dataset/testing.npz')
    A, b = training['A'], training['b']
    A_test, b_test = testing['A'], testing['b']

    print(68 * '*')
    print('Linear Support Vector Machine:')
    print('Logistic loss + Ridge regularizer')
    print('dataset size : {} x {}'.format(A.shape[0], A.shape[1]))
    print(68 * '*')

    print('Numerical solution process is started:')

    # Set parameters and solve numerically with methods from algorithms.py.
    n, p = A.shape
    lmbd = 0.1
    parameter = {}
    parameter['Lips'] = (1 / 2) * LA.norm(A, 'fro')**2 + lmbd
    # parameter['Lips'] = (1 / 2) * LA.norm(A, 'fro')**2 + lmbd, tighter estimate
    parameter['strcnvx'] = lmbd
    parameter['x0'] = np.zeros(p)
    parameter['Lmax'] = 0
    for i in range(n):
        parameter['Lmax'] = np.maximum(parameter['Lmax'],
                                       LA.norm(A[i], 2) * LA.norm(A[i], 2))
    parameter['Lmax'] += lmbd

    parameter['iter_print'] = 100

    fx, gradf, gradfsto = Oracles(b, A, lmbd)
    prox_fx, prox_gradf, prox_gradfsto = Oracles(b, A, 0)
    x, info, error = {}, {}, {}

    # first-order methods
    parameter['maxit'] = 4000
    x['GD'], info['GD'] = GD(fx, gradf, parameter)
    error['GD'] = compute_error(A_test, b_test, x['GD'])
    print('Error w.r.t 0-1 loss: {}'.format(error['GD']))

    x['GDstr'], info['GDstr'] = GDstr(fx, gradf, parameter)
    error['GDstr'] = compute_error(A_test, b_test, x['GDstr'])
    print('Error w.r.t 0-1 loss: {}'.format(error['GDstr']))

    # Estimation of the optimal solution.
    xstar = np.array([
        -0.32775484, 0.48385655, 0.46286856, 0.44371049, 0.33732536,
        0.20245555, 0.69386679, 0.36147147, 0.34303099, -0.09288907
    ])

    ## Into Exercise 2.4

    # Compute the theoretical rates (2.4.a,b)
    diff = np.array([np.linalg.norm(x - xstar) for x in info["GD"]["x"]])
    diffstr = np.array([np.linalg.norm(x - xstar) for x in info["GDstr"]["x"]])

    diffinit = diff[0]  #\|x^0-x^*\|
    diffinit_str = diffstr[0]  #\|x^0-x^*\|
    L = parameter['Lips']
    mu = parameter['strcnvx']
    rate_GD = np.log((L - mu) / (L + mu)) / 2
    rate_GDmu = np.log((L - mu) / (L + mu))
    theoretical_rate_gd = lambda k: ((L - mu) / (L + mu))**(k / 2)
    theoretical_rate_gd_str = lambda k: ((L - mu) / (L + mu))**k

    iter_ = np.arange(len(diff))

    # The plotting function could be plt.semilogy, plt.semilogy, plt.semilogx or plt.plot
    fig, (ax1, ax2) = plt.subplots(1, 2)
    ax1.semilogy(iter_,
                 theoretical_rate_gd(iter_) * diffinit,
                 "k--",
                 label="Theoretical rate GD")
    ax1.semilogy(iter_, diff, "k-", label="Empirical rate GD")
    ax1.legend()
    ax1.grid()
    ax2.semilogy(iter_,
                 theoretical_rate_gd_str(iter_) * diffinit_str,
                 "r--",
                 label="Theoretical rate GD_str")
    ax2.semilogy(iter_, diffstr, "r-", label="Empirical rate GD_str")
    ax2.legend()
    ax2.grid()
    plt.savefig("figs/fig_ex2_4_cconvergence_rates.pdf")

    ## Linear regression (2.4.c)
    # Linear regression: fit the coefficients (slope,intercept) s.t. Y=slope*iter_+intercept
    iter_ = iter_.reshape(iter_.shape[0], 1)
    Y = np.log(diff)
    reg = sklearn.linear_model.LinearRegression().fit(iter_, Y)

    slope_GD = reg.coef_[0]
    intercept_GD = reg.intercept_

    Ymu = np.log(diffstr)
    reg = sklearn.linear_model.LinearRegression().fit(iter_, Ymu)

    slope_GDmu = reg.coef_[0]
    intercept_GDmu = reg.intercept_

    ## Transform back the slope and intercept into the coefficients (a,b)
    a_GD = np.exp(slope_GD)
    b_GD = np.exp(intercept_GD)
    a_GDmu = np.exp(slope_GDmu)
    b_GDmu = np.exp(intercept_GDmu)

    print("GD results: a_GD={:.6f}, b_GD={:.6f}".format(a_GD, b_GD))
    print("GD theoretical rate: a_GD={:.6f}, b_GD={:.6f}".format(
        rate_GD, diffinit))
    print("GD rate diff ={}".format(a_GD - rate_GD))
    print("GD_str results: a_GD={:.6f}, b_GD={:.6f}".format(a_GDmu, b_GDmu))
    print("GD_str theoretical rate: a_GD={:.6f}, b_GD={}".format(
        rate_GDmu, diffinit_str))
    print("GD_str rate diff ={:.6f}".format(a_GDmu - rate_GDmu))

    ## Linear regression_plot (2.4.d)
    fig, (ax1, ax2) = plt.subplots(1, 2)
    ax1.plot(iter_,
             np.log(theoretical_rate_gd(iter_)) + np.log(diffinit),
             "k--",
             label="Theoretical rate GD")
    ax1.plot(iter_, np.log(diff), "k-", label="Empirical rate GD")
    ax1.plot(iter_,
             np.dot(slope_GD, iter_) + intercept_GD,
             "kx-.",
             markevery=500,
             label="Linear regression GD")
    ax1.legend()
    ax1.grid()
    ax2.plot(iter_, np.log(diffstr), "r-", label="Empirical rate GD_str")
    ax2.plot(iter_,
             np.log(theoretical_rate_gd_str(iter_)) + np.log(diffinit_str),
             "r--",
             label="Theoretical rate GD_str")

    ax2.plot(iter_,
             np.dot(slope_GDmu, iter_) + intercept_GDmu,
             "rx-.",
             markevery=500,
             label="Linear regression GD_str")
    ax2.legend()
    ax2.grid()
    plt.savefig("figs/fig_ex2_4_convergence_rates_full.pdf")
    plt.show()


if __name__ == "__main__":
    main()
