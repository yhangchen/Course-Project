import numpy as np


def print_start_message(method_name):
    print('\n\n\n---------- Optimization with {:s} started.\n\n'.format(
        method_name))


def print_end_message(method_name, time_spent):
    print('\n---------- Training over - {:s}. Took {:d} seconds. \n\n'.format(
        method_name, np.math.ceil(time_spent)))


def print_progress(i, maxit, val_F, val_f, val_g):
    print('Iter = {:d}/{:d}, F(X) = {:f}, f(X) = {:f}, g(X) = {:f}'.format(
        i, maxit, val_F, val_f, val_g))
