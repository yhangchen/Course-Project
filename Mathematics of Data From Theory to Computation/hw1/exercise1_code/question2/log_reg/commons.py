import numpy as np


#  compute f
def f(A, b, x):
    fx = np.zeros(A.shape[0])
    for i in range(A.shape[0]):
        ci = b[i] * np.dot(A[i], x)
        fx[i] = np.log(1 + np.exp(-ci))
    fx = fx.mean()
    return fx


# compute gradient
def stogradfx(A, b, x, i):
    di = np.exp(-np.dot(b[i], np.dot(A[i], x)))
    gradfx = (-b[i] * A[i] * di / (1 + di)).T
    return gradfx


def gradfx(A, b, x):
    n = A.shape[0]
    gradfx = np.zeros(x.shape)
    for i in range(n):
        gradfx += stogradfx(A, b, x, i)
    gradfx = gradfx / n

    return gradfx


def Oracles(b, A, lbd):
    """
    FIRST ORDER ORACLE
    Takes inputs b, A, lbd and returns two anonymous functions, one for
    the objective evaluation and the other for the gradient.
    fx(x) computes the objective (l-2 regularized) of input x
    gradf(x) computes the gradient (l-2 regularized) of input x
    gradfsto(x,i) computes the stochastic gradient (l-2 regularized) of input x at index i
    """
    n, p = A.shape
    fx = lambda x: 0.5 * lbd * np.linalg.norm(x, 2)**2 + f(A, b, x)
    gradf = lambda x: lbd * x + gradfx(A, b, x)
    gradfsto = lambda x, i: lbd * x + stogradfx(A, b, x, i)
    return fx, gradf, gradfsto


def compute_error(A_test, b_test, x):
    n, err = A_test.shape[0], 0
    for i in range(n):
        if np.dot(b_test[i], np.dot(A_test[i], x)) <= 0:
            err += 1
    err = err / float(n)
    return err
