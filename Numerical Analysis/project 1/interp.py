import matplotlib.pyplot as plt
from math import cos, pi
plt.style.use("ggplot")
def interp(x, y, method, m0=None):
    """
    This program implememts five integration methods.

    Args:
        :param x: points on the x-axis
        :param y: corresponding values
        :param m0: corresponding differential values
        :param method: specified method

    Returns:
        :interp_fun: Interpolation function.

    """
    assert len(x) == len(y)

    n = len(x)

    if method == 'Hermite':
        assert len(m0) == len(x)

    def interp_fun(t):
        if method == 'Newton':
            diff = [[0 for i in range(n)] for j in range(n)]

            # diff stores the divided difference.
            for i in range(0, n):
                diff[i][i] = y[i]
            for k in range(1, n):
                for i in range(0, n-k):
                    diff[i][i+k] = (diff[i+1][i+k]-diff[i][i+k-1]) / (x[i+k]-x[i])
                    
            ans = diff[0][n-1]
            for i in range(n-2, -1, -1):
                ans = ans*(t-x[i]) + diff[0][i]
            return ans
        
        elif method == 'Lagrange':
            ans = 0
            for i in range(n):
                tmp = y[i]
                for j in range(n):
                    if j != i:
                        tmp *= (t-x[j]) / (x[i]-x[j])
                ans += tmp
            return ans
        
        elif method == 'Linear':
            for i in range(0, n):
                if x[i] <= t <= x[i+1]:
                    ans = (y[i+1]*(t-x[i]) - y[i]*(t-x[i+1])) / (x[i+1]-x[i])
                    return ans
        
        elif method == 'Hermite':
            # m = [-2*t/(1+t**2)**2 for t in x]
            for i in range(n):
                if x[i] <= t <= x[i+1]:
                    delta = x[i+1]-x[i]
                    left = t - x[i]
                    right = t - x[i+1]
                    ans = y[i] * (1+2*left/delta) * (right/delta)**2
                    ans += y[i+1] * (1-2*right/delta) * (left/delta)**2 
                    ans += m0[i] * left * (right/delta)**2
                    ans += m0[i+1] * right * (left/delta)**2
                    return ans
        
        elif method == 'NaturalSpline':
            def Tridiag(sub, dig, upper, b):
                """
                Solve the tridiagnoal system Ax = b.

                Args:
                    :param sub: the lower diagnoal of A
                    :param dig: the main diagnoal of A
                    :param upper: the upper diganoal of A
                    :param b: right terms
                
                Returns:
                    :x: x = A^{-1}b
                """

                m = len(dig)

                assert len(sub) == m-1 and len(upper) == m-1 and len(b) == m

                sub = [0] + sub
                for i in range(1, m):
                    dig[i] -= sub[i]*upper[i-1]/dig[i-1]
                    b[i] -= sub[i]*b[i-1]/dig[i-1]

                x = [0] * m
                x[m-1] = b[m-1]/dig[m-1]
                for i in range(m-2, -1, -1):
                    x[i] = (b[i]-upper[i]*x[i+1])/dig[i]
                return x

            delta = [x[1]-x[0]]
            lamb = [1]
            mu = [3*delta[0]*(y[1]-y[0])]
            
            # Construct the equation matrix
            for i in range(1, n-1):
                delta.append(x[i+1]-x[i])
                lamb.append(delta[i]/(delta[i-1]+delta[i]))
                mu.append(3*(1-lamb[i])/delta[i-1]*(y[i]-y[i-1])+3*lamb[i]/delta[i]*(y[i+1]-y[i]))
            mu.append(3*delta[n-2]*(y[n-1]-y[n-2]))
            m_lamb = [1-t for t in lamb[1:n-1]]
            m_lamb.append(1)

            m = Tridiag(m_lamb, [2] * n, lamb, mu)

            for i in range(n):
                if x[i] <= t <= x[i+1]:
                    delta = x[i+1]-x[i]
                    left = t - x[i]
                    right = t - x[i+1]
                    ans = y[i] * (1+2*left/delta) * (right/delta)**2
                    ans += y[i+1] * (1-2*right/delta) * (left/delta)**2 
                    ans += m[i] * left * (right/delta)**2
                    ans += m[i+1] * right * (left/delta)**2
                    return ans
        else:
            raise ValueError
    return interp_fun



def test(method):
    R = lambda x: 1/(1+x**2)
    if method == 'Lagrange':
        x = [5*cos((2*i+1)*pi/42) for i in range(21)]
    else:
        x = range(-5, 6)
    
    y = [R(t) for t in x]

    if method == 'Hermite':
        m = [-2*t/(1+t**2)**2 for t in x]
    else:
        m = None

    interp_fun = interp(x, y, method, m)

    X = [i/100-5 for i in range(1000)]
    
    plt.grid()
    plt.plot(X, [R(t) for t in X], label="Runge", linewidth=2)
    plt.plot(X, [interp_fun(t) for t in X], label=method)
    plt.legend(loc='upper right')
    plt.xlim(-5.5, 5.5)
    fig = plt.gcf()
    fig.savefig("%s.jpg" % method)
    plt.show()

    plt.grid()
    plt.plot(X, [interp_fun(t)-R(t) for t in X], label='error')
    plt.legend(loc='upper right')
    fig = plt.gcf()
    fig.savefig("%s_error.jpg" % method)
    plt.show()

if __name__ == "__main__":
    methods = ['Newton', 'Lagrange', 'Linear', 'Hermite', 'NaturalSpline']
    for method in methods:
        test(method)
