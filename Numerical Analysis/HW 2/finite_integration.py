from math import sqrt, exp, pi, sin
import numpy as np
import matplotlib.pyplot as plt
plt.style.use("ggplot")

def midpoint(f, L, n):
    ans = 0
    for i in range(n):
        ans += f((i+0.5)*L/n)
    return ans * L /n

def simpson(f, L, n):
    h = L / n
    ans = 2/3*f(h/2)+f(h)/6
    for i in range(1, n):
        ans += f(i*h)/6+2/3*f(i*h+h/2)+f(i*h+h)/6
    return ans * h


def composite(f, L, n, method):
    """
    docstring here

    Input:
        :param f: function on [0, +\infty)
        :param L: length of the integration interval
        :string: method: 'midpoint','trapezoidal' and 'simpson'.
    Return:
        :ans: the value of \int_0^{L} f(x) dx
    """
    h = L/n; eps = 1e-10


    if method == 'midpoint':
        T0 = midpoint(f, L, n)
        T1 = midpoint(f, L, 2*n)
        count = 0
        while abs(T0-T1) > eps:
            count += 1
            T0 = T1
            n *= 2
            T1 = midpoint(f, L, 2*n)
        return T1, count

    elif method == 'trapezoidal':
        T0 = 0; T1 = 0
        for i in range(1, n):
            T0 += f(i*h)
            T1 += f((i-0.5)*h)
        T0 += f(L)/2
        T1 += f(L-0.5*h) + T0
        T0 *= h
        T1 *= h/2
        count = 0
        while abs(T0-T1) > eps:
            T0 = T1
            T1 /= 2
            n *= 2
            h /= 2
            count += 1
            for i in range(n):
                T1 += f((i+0.5)*h)*h/2
        return T1, count

    elif method == 'simpson':
        T0 = simpson(f, L, n)
        T1 = simpson(f, L, 2*n)
        count = 0
        while abs(T0-T1) > eps:
            count += 1
            T0 = T1
            n *= 2
            T1 = midpoint(f, L, 2*n)
        return T1, count    
    else:
        raise ValueError('incorrect type of quadrature')



def Romberg(f, L, n):
    """
    docstring here

    Input:
        :param f: function on [0, +\infty)
        :param L: length of the integration interval
        :param n: maximum number of grid point is 2^n
    Return:
        :ans: the value of \int_0^{L} f(x) dx
    """
    romberg = np.zeros((n,n))
    h = np.divide(L, np.exp2(range(n)))
    romberg[0,0] = L / 2 * f(L)
    for j in range(1,n):
        tmp = 0
        for i in range(2**(j-1)+1):
            tmp += f((2 * i + 1) * h[j])
        romberg[j,0] = romberg[j-1,0] / 2 + tmp * h[j]
        for k in range(1, j+1):
            romberg[j, k] = (4**k * romberg[j,k-1] - romberg[j-1,k-1]) / (4**k - 1)
    return romberg

def Legendre(t, n):
    """
    docstring here
    
    Input:
        :param t: in [-1,1]
        :param n: degree of the Laguerre polynomial

    Return:
        :ans: list, ans[0] = L_n(t); ans[1] = L'n(t)
    """
    
    if n == 0:
        return 1
    elif n == 1:
        return t
    else:
        tmp = [1 for i in range(n+1)]
        tmp[1] = t
        for i in range(2, n+1):
            tmp[i] = (tmp[i-1]*(2*i-1)*t - tmp[i-2]*(i-1)) / i
        return tmp[-1]

def Legendre_coeff(x, n):
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
        l = Legendre(x[i], n+1)
        A[i] = 2*(1-x[i]**2)/((n+1)*l)**2
    return A

def Legendre_zeros(n):
    """
    docstring here

    Input:
        :param n: degree of the Laguerre polynomial
    
    Return:
        :roots: a list containing all the roots of the polynomial
    """
    c = [0 for _ in range(n)] + [1]
    return np.polynomial.legendre.legroots(c)

def Gauss_Legendre(f, L, n):
    """
    docstring here

    Input:
        :param f: function on [0, +\infty)
    
    Return:
        :ans: the value of \int_0^{+\infty} e^{-x}f(x) dx
    """
    x = Legendre_zeros(n)
    A = Legendre_coeff(x, n)
    ans = 0
    for i in range(n):
        ans += f((x[i]+1)*L/2)*A[i]
    return ans*L/2

def newton_cotes_coeff(rn):
    '''
    Parameters
    ----------
    rn : int
        The integer order for equally-spaced data or the relative positions of
        the samples with the first sample at 0 and the last at L, where N+1 is
        the length of `rn`.  N is the order of the Newton-Cotes integration.
    Returns
    -------
    an : ndarray
        1-D array of weights to apply to the function at the provided sample
        positions.
    modified from 
    https://github.com/scipy/scipy/blob/v0.18.0/scipy/integrate/quadrature.py#L833-L842
    '''
    if (rn[0] != 0):
        raise ValueError("The sample positions must start at 0"
                         " and end at N")
    L = rn[-1]
    N = len(rn)-1
    yi = rn / float(L)
    ti = 2 * yi - 1
    nvec = np.arange(N+1)
    C = ti ** nvec[:, np.newaxis]
    Cinv = np.linalg.inv(C)
    # improve precision of result
    for i in range(2):
        Cinv = 2*Cinv - Cinv.dot(C).dot(Cinv)
    vec = 2.0 / (nvec[::2]+1)
    ai = Cinv[:, ::2].dot(vec) * (L / 2.)
    return ai


def newton_cotes(f, L, n):
    """
    docstring here

    Input:
        :param f: function on [0, +\infty)
    
    Return:
        :ans: the value of \int_0^{+\infty} e^{-x}f(x) dx
    """
    rn = np.arange(n+1)*L/n
    A = newton_cotes_coeff(rn)
    ans = 0
    for i in range(n):
        ans += f(float(rn[i]))*A[i]
    return ans
def test(L, n, method):
    def f(x):
        if x < 1e-10:
            return 0
        else:
            return x**3*exp(-x)/(1-exp(-x))
    if method == 'midpoint' or method == 'trapezoidal' or method == 'simpson':
        est, count = composite(f, L, n, method)
        err = abs(est - pi**4/15)
        return err, count
    elif method == 'romberg':
        r = Romberg(f, L, n)
        err = np.abs(np.diag(r) - pi**4/15)
        count = n
        plt.semilogy(range(n), err, marker='x')
        plt.xlabel(r'$\log_2(n)$')
        plt.ylabel('error')
        plt.xticks(range(20))
        plt.title('Romberg')
        import tikzplotlib
        tikzplotlib.save("romberg.tex")
        plt.show()
        return err
    elif method == 'legendre':
        r = Gauss_Legendre(f, L, n)
        err = np.abs(r - pi**4/15)
        return err
    elif method == 'cotes':
        r = newton_cotes(f, L, n)
        err = np.abs(r - pi**4/15)
        return err
    else:
        raise ValueError('wrong method')
    
if __name__ == "__main__":
    # romberg
    test(50,20,'romberg')
    # composite
    for method in ['midpoint', 'trapezoidal', 'simpson']:
        print(method)
        for L in [10,20,40,80,160]:
            print(test(L, 10, method))
    # newton cotes
    for i in [5,10,20,30,40,50]:
        print(test(50,i,'cotes'))

    # Gauss legendre
    err_legendre = []
    for n in [10,15,20,25,30,35,40,45,50,55,60,65,70,75,80]:
        err_legendre.append(test(50, n, 'legendre'))
    plt.semilogy([10,15,20,25,30,35,40,45,50,55,60,65,70,75,80], err_legendre, marker='x')
    plt.xlabel(r'$n$')
    plt.ylabel('error')
    plt.title(r'Legendre, $L=50$')
    import tikzplotlib
    tikzplotlib.save("legendre.tex")
    plt.show()
    err_legendre = []
    for n in [10,15,20,25,30,35,40,45,50,55,60,65,70,75,80]:
        err_legendre.append(test(200, n, 'legendre'))
    plt.semilogy([10,15,20,25,30,35,40,45,50,55,60,65,70,75,80], err_legendre, marker='x')
    plt.xlabel(r'$n$')
    plt.ylabel('error')
    plt.title(r'Legendre, $L = 200$')
    import tikzplotlib
    tikzplotlib.save("legendre2.tex")
    plt.show()