from math import sqrt, exp, pi, factorial
import numpy
import matplotlib.pyplot as plt
plt.style.use("ggplot")

def Laguerre(t, n):
    """
    docstring here
    
    Input:
        :param t: in [0,\infty)
        :param n: degree of the Laguerre polynomial

    Return:
        :ans: list, ans[0] = L_n(t); ans[1] = L'n(t)
    """
    
    if n == 0:
        return [1,0]
    elif n == 1:
        return [1-t,-1]
    else:
        tmp = [1 for i in range(n+1)]
        tmp[1] = 1-t
        for i in range(2, n+1):
            tmp[i] = (tmp[i-1]*(2*i-1-t) - tmp[i-2]*(i-1)) / i
        diff = n * (tmp[-1] - tmp[-2]) / t
        return [tmp[-1], diff]

def Laguerre_coeff(x, n):
    """
    docstring here

    Input:
        :param x: in [0,\infty)
        :param n: degree of the Laguerre polynomial

    Return:
        :A: a list containing the parameters.
    """ 
    A = [0 for _ in range(n)]
    for i in range(n):
        l = Laguerre(x[i], n+1)[0]
        A[i] = x[i]/((n+1)*l)**2
    return A

def Laguerre_zeros(n):
    """
    docstring here

    Input:
        :param n: degree of the Laguerre polynomial
    
    Return:
        :roots: a list containing all the roots of the polynomial
    """
    c = [0 for _ in range(n)] + [1]
    return numpy.polynomial.laguerre.lagroots(c)

def Gauss_Laguerre(f, n):
    """
    docstring here

    Input:
        :param f: function on [0, +\infty)
    
    Return:
        :ans: the value of \int_0^{+\infty} e^{-x}f(x) dx
    """
    x = Laguerre_zeros(n)
    A = Laguerre_coeff(x, n)
    ans = 0
    for i in range(n):
        ans += f(x[i])*A[i]
    return ans



def test(n, method):
    if method == 'Laguerre':
        f = lambda x: x**3/(1-exp(-x))
        err = Gauss_Laguerre(f, n) - pi**4/15
        return abs(err)


if __name__ == "__main__":
    '''
    proceed until overflow
    '''
    n = [i+1 for i in range(200)]
    err = [test(n[i], 'Laguerre') for i in range(200)]
    plt.loglog(n, err, marker='x')
    plt.xlabel(r'$n$')
    plt.ylabel('error')
    plt.title('Gauss-Laguerre')
    import tikzplotlib
    tikzplotlib.save("gauss.tex")
    plt.show()
