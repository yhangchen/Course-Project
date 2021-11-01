from itertools import cycle
from collections import defaultdict
from typing import ItemsView

import numpy as np
from numpy import linalg as LA
from log_reg.commons import Oracles, compute_error, gradfx
from log_reg.algorithms import GD, GDstr, AGD, AGDstr, AGDR, AdaGrad, SGD, SAG, SVR
# from log_reg.algorithms import LSGD, LSAGD, LSAGDR, ADAM
from log_reg.algorithms import SubG, fista, ista, prox_sg
from log_reg.operators import l1_prox, l2_prox, norm1, norm2sq
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors


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

    # Choose the solvers you want to call
    methods: defaultdict[str, bool] = defaultdict(lambda: False)
    methods['GD'] = True
    methods['GDstr'] = True
    methods['AGD'] = True
    methods['AGDstr'] = True
    methods['LSGD'] = False
    methods['LSAGD'] = False
    methods['AGDR'] = True
    methods['LSAGDR'] = False
    methods['AdaGrad'] = True
    methods['ADAM'] = False
    methods['SGD'] = True
    methods['SAG'] = True
    methods['SVR'] = True

    methods['SubG'] = True
    methods['ISTA'] = True
    methods['FISTA'] = True
    methods['FISTAR'] = True
    methods['PROXSG'] = True

    # Choose what regularizers to use for the prox
    prox_regularizers = ['L1', 'L2']

    print('Numerical solution process is started:')

    # Set parameters and solve numerically with methods from algorithms.py.
    n, p = A.shape
    lmbd = 0.1
    f_star = 0.278635
    parameter = {}
    parameter['Lips'] = (1 / 2) * LA.norm(A, 'fro')**2 + lmbd
    parameter['strcnvx'] = lmbd
    parameter['x0'] = np.zeros(p)
    parameter['Lmax'] = 0
    for i in range(n):
        parameter['Lmax'] = np.maximum(parameter['Lmax'],
                                       LA.norm(A[i], 2) * LA.norm(A[i], 2))
    parameter['Lmax'] += lmbd

    # prox config
    f_star_l1 = 0.675361
    f_star_l2 = 0.278635
    lmbd_l1 = 0.3
    lmbd_l2 = 0.1
    parameter['iter_print'] = 100
    parameter['prox_Lips'] = (1 / 2) * LA.norm(A, 'fro')**2
    parameter['stoch_rate_regime'] = lambda k: 1 / (parameter['prox_Lips'] + k
                                                    **(0.55))

    # Stochastic methods
    parameter['no0functions'] = n

    fx, gradf, gradfsto = Oracles(b, A, lmbd)
    prox_fx, prox_gradf, prox_gradfsto = Oracles(b, A, 0)
    x, info, error = {}, {}, {}

    # first-order methods
    parameter['maxit'] = 4000
    if methods['GD']:
        x['GD'], info['GD'] = GD(fx, gradf, parameter)
        error['GD'] = compute_error(A_test, b_test, x['GD'])
        print('Error w.r.t 0-1 loss: {}'.format(error['GD']))

    parameter['maxit'] = 4000
    if methods['GDstr']:
        x['GDstr'], info['GDstr'] = GDstr(fx, gradf, parameter)
        error['GDstr'] = compute_error(A_test, b_test, x['GDstr'])
        print('Error w.r.t 0-1 loss: {}'.format(error['GDstr']))

    parameter['maxit'] = 5000
    if methods['AGD']:
        x['AGD'], info['AGD'] = AGD(fx, gradf, parameter)
        error['AGD'] = compute_error(A_test, b_test, x['AGD'])
        print('Error w.r.t 0-1 loss: {}'.format(error['AGD']))

    parameter['maxit'] = 5000
    if methods['AGDstr']:
        x['AGDstr'], info['AGDstr'] = AGDstr(fx, gradf, parameter)
        error['AGDstr'] = compute_error(A_test, b_test, x['AGDstr'])
        print('Error w.r.t 0-1 loss: {}'.format(error['AGDstr']))

    parameter['maxit'] = 300
    if methods['AGDR']:
        x['AGDR'], info['AGDR'] = AGDR(fx, gradf, parameter)
        error['AGDR'] = compute_error(A_test, b_test, x['AGDR'])
        print('Error w.r.t 0-1 loss: {}'.format(error['AGDR']))

    parameter['maxit'] = 400
    if methods['LSGD']:
        x['LSGD'], info['LSGD'] = LSGD(fx, gradf, parameter)
        error['LSGD'] = compute_error(A_test, b_test, x['LSGD'])
        print('Error w.r.t 0-1 loss: {}'.format(error['LSGD']))

    if methods['LSAGD']:
        x['LSAGD'], info['LSAGD'] = LSAGD(fx, gradf, parameter)
        error['LSAGD'] = compute_error(A_test, b_test, x['LSAGD'])
        print('Error w.r.t 0-1 loss: {}'.format(error['LSAGD']))

    parameter['maxit'] = 30
    if methods['LSAGDR']:
        x['LSAGDR'], info['LSAGDR'] = LSAGDR(fx, gradf, parameter)
        error['LSAGDR'] = compute_error(A_test, b_test, x['LSAGDR'])
        print('Error w.r.t 0-1 loss: {}'.format(error['LSAGDR']))

    parameter['maxit'] = 4000
    if methods['AdaGrad']:
        x['AdaGrad'], info['AdaGrad'] = AdaGrad(fx, gradf, parameter)
        error['AdaGrad'] = compute_error(A_test, b_test, x['AdaGrad'])
        print('Error w.r.t 0-1 loss: {}'.format(error['AdaGrad']))

    parameter['maxit'] = 4000
    if methods['ADAM']:
        x['ADAM'], info['ADAM'] = ADAM(fx, gradf, parameter)
        error['ADAM'] = compute_error(A_test, b_test, x['ADAM'])
        print('Error w.r.t 0-1 loss: {}'.format(error['ADAM']))

    # stochastic methods
    parameter['maxit'] = 5 * n
    if methods['SGD']:
        x['SGD'], info['SGD'] = SGD(fx, gradfsto, parameter)
        error['SGD'] = compute_error(A_test, b_test, x['SGD'])
        print('Error w.r.t 0-1 loss: {}'.format(error['SGD']))

    if methods['SAG']:
        x['SAG'], info['SAG'] = SAG(fx, gradfsto, parameter)
        error['SAG'] = compute_error(A_test, b_test, x['SAG'])
        print('Error w.r.t 0-1 loss: {}'.format(error['SAG']))

    parameter['maxit'] = int(.2 * n)
    if methods['SVR']:
        x['SVR'], info['SVR'] = SVR(fx, gradf, gradfsto, parameter)
        error['SVR'] = compute_error(A_test, b_test, x['SVR'])
        print('Error w.r.t 0-1 loss: {}'.format(error['SVR']))

    # plot figures for non-prox methods
    print('Numerical solution process of non-prox methods is completed.')
    det_methods = [
        'GD', 'GDstr', 'AGD', 'AGDstr', 'LSGD', 'LSAGD', 'AGDR', 'LSAGDR',
        'AdaGrad', 'ADAM'
    ]
    stoc_methods = ['GD', 'SGD', 'SAG', 'SVR']
    plot_methods(x, info, "figs/fig_ex1.pdf", f_star, n, methods, det_methods,
                 stoc_methods)

    # Prox methods
    for reg in prox_regularizers:
        if reg == 'L1':
            gx = norm1
            prox_fc = l1_prox
            prox_f_star = f_star_l1
            parameter['lambda'] = lmbd_l1
        elif reg == 'L2':
            gx = norm2sq
            prox_fc = l2_prox
            prox_f_star = f_star_l2
            parameter['lambda'] = lmbd_l2
        else:
            raise ValueError("Prox regularizer not supported")

        parameter['maxit'] = 5000
        if methods['SubG'] and reg == 'L1':
            x['SubG'], info['SubG'] = SubG(prox_fx, gx, prox_gradf, parameter)
            error['SubG'] = compute_error(A_test, b_test, x['SubG'])
            print('Error w.r.t 0-1 loss: {}'.format(error['SubG']))

        if methods['ISTA']:
            x['ISTA'], info['ISTA'] = ista(prox_fx, gx, prox_gradf, prox_fc,
                                           parameter)
            error['ISTA'] = compute_error(A_test, b_test, x['ISTA'])
            print('Error w.r.t 0-1 loss: {}'.format(error['ISTA']))

        parameter['maxit'] = 2000
        if methods['FISTA']:
            parameter['restart_fista'] = False
            x['FISTA'], info['FISTA'] = fista(prox_fx, gx, prox_gradf, prox_fc,
                                              parameter)
            error['FISTA'] = compute_error(A_test, b_test, x['FISTA'])
            print('Error w.r.t 0-1 loss: {}'.format(error['FISTA']))

        if methods['FISTAR']:
            parameter['restart_fista'] = True
            x['FISTAR'], info['FISTAR'] = fista(prox_fx, gx, prox_gradf,
                                                prox_fc, parameter)
            error['FISTAR'] = compute_error(A_test, b_test, x['FISTAR'])
            print('Error w.r.t 0-1 loss: {}'.format(error['FISTAR']))

        parameter['maxit'] = 5000
        if methods['PROXSG']:
            x['PROXSG'], info['PROXSG'] = prox_sg(prox_fx, gx, prox_gradfsto,
                                                  prox_fc, parameter)
            error['PROXSG'] = compute_error(A_test, b_test, x['PROXSG'])
            print('Error w.r.t 0-1 loss: {}'.format(error['PROXSG']))

        # Plot prox methods
        print(
            f'Numerical solution process of {reg} prox methods is completed.')
        if reg == 'L1':
            det_methods = ['SubG', 'ISTA', 'FISTA', 'FISTAR']
        else:
            det_methods = ['ISTA', 'FISTA', 'FISTAR']
        stoc_methods = ['SubG', 'ISTA', 'PROXSG']
        plot_methods(x, info, f'figs/fig_ex1_prox_{reg}.pdf', prox_f_star, n,
                     methods, det_methods, stoc_methods)
        for key, value in error.items():
            print(key, ':', value)


def plot_methods(x,
                 info,
                 path,
                 f_star,
                 n,
                 methods,
                 det_methods=[],
                 stoc_methods=[]):
    # Only plot if some of the selected methods have been run
    if any([methods[name] for name in det_methods + stoc_methods]):

        colors = cycle([(0, 0, 1), (0, 0.5, 0), (1, 0, 0), (0, 0.75, 0.75),
                        (0.75, 0, 0.75), (0.75, 0.75, 0), (0, 0, 0),
                        (0.5, 1, 0.5), (0.5, 0, 0.5), (0.75, 0.25, 0.25)])

        ax1 = plt.subplot(1, 2, 1)

        for key in x.keys():
            if key in det_methods:
                ax1.plot(np.array(range(info[key]['iter'])),
                         info[key]['fx'] - f_star,
                         color=next(colors),
                         lw=2,
                         label=key)
        ax1.legend()
        #ax1.set_ylim(1e-9, 1e0)
        ax1.set_xlabel('#iterations')
        ax1.set_ylabel(r'$f(\mathbf{x}^k) - f^\star$')
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.grid()

        ax2 = plt.subplot(1, 2, 2)

        for key in x.keys():
            if key in stoc_methods:
                if key in ['GD', 'ISTA']:
                    ax2.plot(np.array(range(info[key]['iter'])),
                             np.maximum(info[key]['fx'] - f_star, 0),
                             lw=2,
                             color=next(colors),
                             label=key,
                             marker='o')
                else:
                    ax2.plot(np.array(range(info[key]['iter'])) / float(n),
                             np.maximum(info[key]['fx'] - f_star, 0),
                             lw=2,
                             color=next(colors),
                             label=key)
        ax2.set_xlim(0, 5)
        #ax2.set_ylim(1e-4, 1e0)
        ax2.legend()
        ax2.set_xlabel('#epochs')
        ax2.set_ylabel(r'$f(\mathbf{x}^k) - f^\star$')
        ax2.set_yscale('log')
        ax2.grid()

        plt.tight_layout()
        plt.savefig(path)
        plt.clf()


if __name__ == "__main__":
    main()
