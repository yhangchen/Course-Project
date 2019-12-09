import numpy as np
from sklearn import linear_model
import matplotlib.pyplot as plt
from scipy import optimize
import time
class Potts(object):
    def __init__(self, args, method = 'MH'):
        # d means the dimension of the Potts model.
        self._q = args['q']
        self._N = args['N']
        self._J = args['J']
        self._h = args['h']
        self._d = args['d']
        self._T = args['T']
        self._coordinates = args['coordinates']
        self._size = args['size']
        self._concrete = np.random.randint(1, self._q + 1, size = self._size)           

    def rand_coordinate(self):
        out = [0 for _ in range(self._d)]
        for i in range(self._d):
            out[i] = np.random.randint(self._N)
        return tuple(out)
    
    def element(self, coordinate):
        return self._concrete[coordinate]
    
    def coordinates(self):
        return self._coordinates
    
    def concrete(self):
        return self._concrete
    
    def flip(self, coordinate, new):
        self._concrete[coordinate] = new

    def neighbour(self, coordinate):
        #output the value in the neighbourhood of coordinate (tuple)
        neighbours = []
        for i in range(self._d):
            p = list(coordinate)
            p[i] = (p[i] - 1) % self._N
            neighbours.append(self._concrete[tuple(p)])
            p[i] = (p[i] + 2) % self._N
            neighbours.append(self._concrete[tuple(p)])
        return neighbours

    def neighbour_coordinate(self, coordinate):
        #output the coordinate of the neighbourhood of coordinate (tuple)
        neighbours = []
        for i in range(self._d):
            p = list(coordinate)
            p[i] = (p[i] - 1) % self._N
            neighbours.append(tuple(p))
            p[i] = (p[i] + 2) % self._N
            neighbours.append(tuple(p))
        return neighbours

    def hamiltonian(self):
        #periodic boundary condition
        concrete = self._concrete
        delta = 0
        for coordinate in self._coordinates:
            neighbour = self.neighbour(coordinate)
            current = concrete[coordinate]
            for element in neighbour:
                delta += int(element==current)
        return -self._J * delta / 2 - self._h * concrete.sum()

    def magnetization(self):
        return self._concrete.sum()
    
    def diff_hamiltonian(self, coordinate, new):
        #flip the coordinate to new, return H(new) - H(old)
        concrete = self._concrete
        old = concrete[coordinate]
        neighbour = self.neighbour(coordinate)
        old_relation = new_relation = 0
        for element in neighbour:
            old_relation += int(element==old)
            new_relation += int(element==new)
        return self._h * (old - new) + self._J * (old_relation - new_relation)

    def correlation_estimate(self, k):
        out = 0
        for coordinate in self._coordinates:
            current = self._concrete[coordinate]
            k_neighbour = list(coordinate)
            for i in range(self._d):
                k_neighbour[i] = (k_neighbour[i] + k) % self._N
                out += current * self._concrete[tuple(k_neighbour)]
                k_neighbour[i] = (k_neighbour[i] - 2*k) % self._N
                out += current * self._concrete[tuple(k_neighbour)]
                k_neighbour[i] = (k_neighbour[i] + k) % self._N
        return out / self._d / 2
    
    def diff_correlation(self, coordinate, k, new):
        k_neighbour = list(coordinate)
        current = self._concrete[coordinate]
        out = 0
        for i in range(self._d):
            k_neighbour[i] = (k_neighbour[i] + k) % self._N
            k_value = self._concrete[tuple(k_neighbour)]
            out += (new - current) * k_value
            k_neighbour[i] = (k_neighbour[i] - 2 * k) % self._N
            k_value = self._concrete[tuple(k_neighbour)]
            out += (new - current) * k_value
            k_neighbour[i] = (k_neighbour[i] + k) % self._N
        return out / self._d / 2


def parameter_specify(q, N, d):
    args = dict()
    args['N'] = N; args['q'] = q; args['d'] = d
    out = [i for i in range(d)]
    coordinates = [0 for _ in range(N**d)]
    for i in range(N**d):
        for j in range(d):
            out[j] = (i//(N**j)) % N
        coordinates[i]= tuple(out)
    args['coordinates'] = coordinates
    args['size'] = tuple(N for _ in range(d))
    return args

    
def Metropolis(args: dict):

    converged = False
    T = args['T']; kB = args['kB']; q = args['q']; N = args['N']; d = args['d']; critical = args['critical']
    k_lst = args['k_lst']; num_k = len(k_lst)
    cor = [0 for _ in range(num_k)]
    beta = 1/(kB*T)
    states = set(i for i in range(1, q+1))

    M =  0
    H =  0
    H_square = 0
    new_average_M =  0
    average_M =  0
    new_average_H =  0
    average_H =  0
    new_average_H_2 =  0
    average_H_2 =  0
    cor =  [0 for _ in range(num_k)]
    average_cor =  np.array([0 for _ in range(num_k)])
    new_average_cor =  np.array([0 for _ in range(num_k)])

    k0 = 0 

    Hamiltonian_lst =  0
    Magnetization_lst =  0
    Correlation_lst = [0 for _ in range(num_k)]

    worker = Potts(args)
    if T < critical: # smaller T requires more warming up steps 
        warming_up = 500000
    else:
        warming_up = 50000
    for _ in range(warming_up):
        # workerropose
        coordinate_flip = worker.rand_coordinate()
        old = worker.element(coordinate_flip)
        new_choices = states - set([old])
        new = np.random.choice(list(new_choices))
        # Decide
        diff = worker.diff_hamiltonian(coordinate_flip, new)
        threshold = min(1, np.exp(-beta*diff))
        r = np.random.random()
        if r <= threshold:
            # flip
            worker.flip(coordinate_flip, new)

    Magnetization_lst = worker.magnetization()
    Hamiltonian_lst = worker.hamiltonian()
    for j in range(num_k):
        Correlation_lst[j] = worker.correlation_estimate(k_lst[j])

    while not converged:
        average_M = new_average_M
        average_H = new_average_H
        average_H_2 = new_average_H_2
        average_cor = new_average_cor.copy()
        M += Magnetization_lst
        h0 = Hamiltonian_lst
        H += h0; H_square += h0**2
        for j in range(num_k):
            cor[j] += Correlation_lst[j]
        k0 += 1
        new_average_M = M / k0; new_average_H = H / k0
        new_average_H_2 = H_square / k0
        for j in range(num_k):
            new_average_cor[j] = cor[j] / k0

        # Propose
        coordinate_flip = worker.rand_coordinate()
        old = worker.element(coordinate_flip)
        new_choices = states - set([old])
        new = np.random.choice(list(new_choices))
        # Decide
        diff = worker.diff_hamiltonian(coordinate_flip, new)
        threshold = min(1, np.exp(-beta*diff))
        r = np.random.random()
        if r <= threshold:
            Magnetization_lst += new - old
            Hamiltonian_lst += diff
            for j in range(num_k):
                Correlation_lst[j] += worker.diff_correlation(coordinate_flip, k_lst[j], new)
            # flip
            worker.flip(coordinate_flip, new)

        err_cor = np.sum(np.abs(new_average_cor - average_cor)) / num_k      

        Err = err_cor + np.abs(new_average_H_2-average_H_2 - (new_average_H**2-average_H**2)) / T**2 + np.abs(new_average_H-average_H) + \
        np.abs(new_average_M - average_M)
        if k0 % 5000 == 0:
            # print('Iter: %s, Error %s' % (k0,Err))
            pass
        elif k0 >= 5*1e3 and Err < 0.001:
            converged = True
        else:
            converged = False
        if k0 >= 1e7:
            converged = True

    out1 = new_average_M / N**d
    out2 = new_average_H / N**d
    Var = new_average_H_2 - new_average_H**2
    out3 = Var / N**d * kB * beta**2
    out4 = new_average_cor / N**d - out1**2 * np.ones(len(new_average_cor))
    return [out1, out2, out3, np.abs(out4)]


def Wolff(args: dict):
    converged = False
    T = args['T']; kB = args['kB']; q = args['q']; N = args['N']; d = args['d']; critical = args['critical']; J = args['J']
    k_lst = args['k_lst']; num_k = len(k_lst)
    cor = [0 for _ in range(num_k)]
    beta = 1/(kB*T)
    states = set(i for i in range(1, q+1))

    M =  0
    H =  0
    H_square = 0
    new_average_M =  0
    average_M =  0
    new_average_H =  0
    average_H =  0
    new_average_H_2 =  0
    average_H_2 =  0
    cor =  [0 for _ in range(num_k)]
    average_cor =  np.array([0 for _ in range(num_k)])
    new_average_cor =  np.array([0 for _ in range(num_k)])

    k0 = 0 

    Hamiltonian_lst =  0
    Magnetization_lst =  0
    Correlation_lst = [0 for _ in range(num_k)]

    worker = Potts(args)

    Magnetization_lst = worker.magnetization()
    Hamiltonian_lst = worker.hamiltonian()
    for j in range(num_k):
        Correlation_lst[j] = worker.correlation_estimate(k_lst[j])
    if T < critical:
        warming_up = 500000
    else:
        warming_up = 10000
    for _ in range(warming_up):
        C_lst = []; new_lst = []
        base = worker.rand_coordinate()
        old = worker.element(base)
        new_lst.append(base)
        checked = set([base])
        while new_lst:
            neonew_lst = []; neighbours = set()
            for current in new_lst:
                neighbourhood = worker.neighbour_coordinate(current)
                for neighbour in neighbourhood:
                    if neighbour not in checked and worker.element(neighbour) == old:
                        r = np.random.random()
                        bond = 1 - np.exp(-beta*J)
                        if r <= bond:
                            neonew_lst.append(neighbour)
                            checked.add(neighbour)
                    neighbours.add(neighbour)
            checked.union(neighbours)
            C_lst.extend(new_lst)
            new_lst = neonew_lst
        C_lst.extend(new_lst)

        new = np.random.choice(list(states-set([old])))

        for i in range(len(C_lst)):
            worker.flip(C_lst[i], new)
    
    while not converged:
        average_M = new_average_M
        average_H = new_average_H
        average_H_2 = new_average_H_2
        average_cor = new_average_cor.copy()

        h0 = Hamiltonian_lst
        H += h0; H_square += h0**2
        M += Magnetization_lst
        for j in range(num_k):
            cor[j] += Correlation_lst[j]

        C_lst = []; new_lst = []
        base = worker.rand_coordinate()
        old = worker.element(base)
        new_lst.append(base)
        checked = set([base])
        while new_lst:
            neonew_lst = []; neighbours = set()
            for current in new_lst:
                neighbourhood = worker.neighbour_coordinate(current)
                for neighbour in neighbourhood:
                    if neighbour not in checked and worker.element(neighbour) == old:
                        r = np.random.random()
                        bond = 1 - np.exp(-beta*J)
                        if r <= bond:
                            neonew_lst.append(neighbour)
                            checked.add(neighbour)
                    neighbours.add(neighbour)
            checked.union(neighbours)
            C_lst.extend(new_lst)
            new_lst = neonew_lst
        C_lst.extend(new_lst)

        new = np.random.choice(list(states-set([old])))
        # for i in range(len(C_lst)):
        #     new[i] = np.random.choice(list(states-set([old])))

        # if len(C_lst) < N**2/4:
        #     for i in range(len(C_lst)):
        #         Hamiltonian_lst += worker.diff_hamiltonian(C_lst[i], new)
        #         Magnetization_lst += new - old
        #         for j in range(num_k):
        #             Correlation_lst[j] += worker.diff_correlation(C_lst[i], k_lst[j], new)
        #         worker.flip(C_lst[i], new)
        # else:
        for i in range(len(C_lst)):
            worker.flip(C_lst[i], new)
        Hamiltonian_lst = worker.hamiltonian()
        Magnetization_lst = worker.magnetization()
        for j in range(num_k):
            Correlation_lst[j] = worker.correlation_estimate(k_lst[j])

        k0 += 1
        new_average_M = M / k0; new_average_H = H / k0
        new_average_H_2 = H_square / k0
        for j in range(num_k):
            new_average_cor[j] = cor[j] / k0
        err_cor = np.sum(np.abs(new_average_cor - average_cor)) / num_k


        Err = err_cor + np.abs(new_average_H_2-average_H_2 - (new_average_H**2-average_H**2)) / T**2 + np.abs(new_average_H-average_H) + \
        np.abs(new_average_M - average_M)
        if k0 % 5000 == 0:
            print('Iter: %s, Error %s' % (k0,Err))
            pass
        if k0 >= 10000 and abs(Err) < 0.001:
            converged = True
        else:
            converged = False
        if k0 >= 100000:
            converged = True


    out1 = new_average_M / N**d
    out2 = new_average_H / N**d
    Var = new_average_H_2 - new_average_H**2
    out3 = Var / N**d * kB * beta**2
    out4 = new_average_cor / N**d - (out1)**2 * np.ones(len(new_average_cor))
    return [out1, out2, out3, np.abs(out4)]





def main_simulation(choice, method):
    J = 1; kB = 1; q = 3; n = 10
    if choice[0:1] == '2':
        d = 2; critical = 1
    elif choice[0:1] == '3':
        d = 3; critical = 1.82
    else:
        raise(ValueError)
    assert method == 'MH' or method == 'Wolff'
    if choice[1:4] == 'Mag':
        if method == 'MH':
            num_h = 15; n = 5
            assert (num_h - 1) % 2 == 0

            if choice[0:1] == '2':
                num_T = 41; T_1 = 0.8; T_2 = 2.0; critical = 1; N = 15
            else:
                num_T = 21; T_1 = 1.5; T_2 = 3.0; critical = 1.82; N = 6

            args = parameter_specify(q, N, d)
            args['J'] = J; args['kB'] = kB; args['k_lst'] = [1]; args['critical'] = critical
            T_lst = np.linspace(T_1, T_2, num_T)
            h_lst = np.linspace(-2,2,num_h)
            Mag_lst = [0 for _ in range(num_T)]
            zero_h = []

            for i in range(num_T):
                T = T_lst[i]
                args['T'] = T
                cache = list()
                print(i)
                for j in range(num_h):
                    args['h'] = h_lst[j]
                    print(j)
                    result = 0
                    ####
                    for i in range(n):
                        ans = Metropolis(args)
                        result += ans[0]
                    result /= n
                    ####
                    cache.append(result)
                    if j == (num_h-1)//2:# h = 0
                        zero_h.append(result)
                Mag_lst[i] = cache
                plt.plot(h_lst, Mag_lst[i], label = 'T = %s' % round(T,2))

            print(Mag_lst)
            print('\n')
            print(zero_h)

            plt.xlabel('h')
            plt.ylabel('Magnetization')
            plt.legend()
            plt.show()

            for i in range(num_T):
                if T_lst[i] <= critical:
                    i += 1
                else:
                    break

            T_lst = T_lst[i+1:]
            zero_h = zero_h[i+1:]
            zero_h = np.array(zero_h)

            log_epsilon = np.log(np.abs(T_lst / critical - np.ones(len(T_lst))))
            log_m = np.log(np.abs(np.array(zero_h)))


            regr = linear_model.LinearRegression()
            regr.fit(log_epsilon.reshape(-1, 1), log_m)
            print('Score: %s' % regr.score(log_epsilon.reshape(-1, 1), log_m))
            print('Coefficient: %s' % regr.coef_)

            plt.subplot(1,2,1)
            plt.plot(log_epsilon, log_m, linewidth = 0, marker = 'o', markersize = 8)
            plt.plot(log_epsilon, regr.predict(log_epsilon.reshape(-1, 1)), color='red', linewidth=2)
            plt.xlabel(r'$\log \epsilon$')
            plt.ylabel(r'$\log m$')
            plt.title('Regression')


            plt.subplot(1,2,2)
            plt.plot(T_lst, zero_h, linewidth = 0, marker = 'o', markersize = 8)
            plt.plot(T_lst, np.exp(regr.predict(log_epsilon.reshape(-1, 1))), color='red', linewidth=2)
            plt.xlabel(r'$T$')
            plt.ylabel(r'$m$')
            plt.title('Prediction')
            plt.show()
        else:
            if choice[0:1] == '2':
                num_T = 21; T_1 = 0.85; T_2 = 1.25; critical = 1; N = 20
            else:
                num_T = 11; T_1 = 1.5; T_2 = 3.0; critical = 1.82; N = 12

            args = parameter_specify(q, N, d)
            args['J'] = J; args['kB'] = kB; args['k_lst'] = [1]; args['critical'] = critical; args['h'] = 0
            T_lst = np.linspace(T_1, T_2, num_T)
            Mag_lst = [0 for _ in range(num_T)]
            zero_h = []

            for i in range(num_T):
                T = T_lst[i]
                print(i)
                args['T'] = T
                ans = Wolff(args)
                result = ans[0]
                zero_h.append(result)
            
            

            log_epsilon = np.log(np.abs(T_lst / critical - np.ones(len(T_lst))))
            log_m = np.log(np.abs(np.array(zero_h)))


            regr = linear_model.LinearRegression()
            regr.fit(log_epsilon.reshape(-1, 1), log_m)
            print('Score: %s' % regr.score(log_epsilon.reshape(-1, 1), log_m))
            print('Coefficient: %s' % regr.coef_)

            plt.subplot(1,2,1)
            plt.plot(log_epsilon, log_m, linewidth = 0, marker = 'o', markersize = 8)
            plt.plot(log_epsilon, regr.predict(log_epsilon.reshape(-1, 1)), color='red', linewidth=2)
            plt.xlabel(r'$\log \epsilon$')
            plt.ylabel(r'$\log m$')
            plt.title('Regression')


            plt.subplot(1,2,2)
            plt.plot(T_lst, zero_h, linewidth = 0, marker = 'o', markersize = 8)
            plt.plot(T_lst, np.exp(regr.predict(log_epsilon.reshape(-1, 1))), color='red', linewidth=2)
            plt.xlabel(r'$T$')
            plt.ylabel(r'$m$')
            plt.title('Prediction')
            plt.show()




    elif choice[1:] == "Hal":
        if method == 'MH':
            n = 10
        else:
            n = 1
        if choice[0:1] == '2':
            num_T = 21; T_lst = np.linspace(0.55, 0.95, num_T); N = 20; critical = 1
        elif choice[0:1] == '3':
            num_T = 15; T_lst = np.linspace(1.2, 3, num_T); N = 8; critical = 1.82
        else:
            raise(ValueError)

        num_k = 1
        args = parameter_specify(q, N, d)
        k_lst = np.array([i for i in range(1, num_k+1)])
        args['J'] = J; args['kB'] = kB; args['n'] = n; args['k_lst'] = k_lst; args['h'] = 0; args['critical'] = critical
        H_lst = [0 for _ in range(num_T)]; H2_lst = [0 for _ in range(num_T)]
        for i in range(num_T):
            args['T'] = T_lst[i]
            H_lst[i] = H2_lst[i] = 0
            print(i)
            for j in range(n):
                if method == 'MH':
                    result = Metropolis(args)
                elif method == 'Wolff':
                    result = Wolff(args)

                H_lst[i] += result[1]
                H2_lst[i] += result[2]
                    
            H_lst[i] /= n
            H2_lst[i] /= n
            print(H2_lst[i])
        print(H_lst)
        print(H2_lst)
        plt.plot(T_lst, H_lst)
        plt.xlabel('T')
        plt.ylabel('Internal Energy')
        plt.show()

        plt.plot(T_lst, H2_lst)
        plt.xlabel('T')
        plt.ylabel('Specific heat')
        plt.show()

        for i in range(num_T):
            if T_lst[i] <= critical:
                i += 1
            else:
                break

        T_lst = T_lst[i+1:]
        H_lst = H_lst[i+1:]
        H2_lst = H2_lst[i+1:]

        epsilon_log = np.log(np.abs(T_lst/critical - np.ones(shape=(len(T_lst)))))
        print(epsilon_log)
        H_lst_log = np.log(np.abs(np.array(H_lst)))
        H2_lst_log = np.log(np.array(H2_lst))

        regr = linear_model.LinearRegression()
        regr.fit(epsilon_log.reshape(-1, 1), H_lst_log)
        print('Score: %s' % regr.score(epsilon_log.reshape(-1, 1), H_lst_log))
        print('Coefficient: %s' % regr.coef_)
        print(epsilon_log)
        print(H_lst_log)
        plt.plot(epsilon_log, H_lst_log, linewidth = 0, marker = 'o', markersize = 8)
        plt.plot(epsilon_log, regr.predict(epsilon_log.reshape(-1, 1)), color='red', linewidth=2)
        plt.xlabel(r'$\log(\epsilon)$')
        plt.ylabel(r'$\log(u)$')
        plt.show()

        regr.fit(epsilon_log.reshape(-1, 1), H2_lst_log)
        print('Score: %s' % regr.score(epsilon_log.reshape(-1, 1), H2_lst_log))
        print('Coefficient: %s' % regr.coef_)
        plt.plot(epsilon_log, H2_lst_log, linewidth = 0, marker = 'o', markersize = 8)
        plt.plot(epsilon_log, regr.predict(epsilon_log.reshape(-1, 1)), color='red', linewidth=2)
        plt.xlabel(r'$log(\epsilon)$')
        plt.ylabel(r'$\log(c)$')
        plt.show()

    elif choice[1:] == 'Cor':
        if method == 'MH':
            n = 10
        else:
            n = 1
        if choice[0:1] =='2':  
            N = 20; num_k = 8; num_T = 10; T_1 = 1.02; T_2 = 1.2
        else:
            N = 12; num_k = 5; num_T = 16; T_1 = 1.85; T_2 = 2.15
        args = parameter_specify(q, N, d)
        k_lst = np.array([i for i in range(1, num_k+1)])
        args['J'] = J; args['kB'] = kB; args['n'] = n; args['k_lst'] = k_lst; args['h'] = 0; args['critical'] = critical
        T_lst = np.linspace(T_1, T_2, num_T)
        Cor_lst = [[0 for _ in range(num_k)] for _ in  range(num_T)]
        for i in range(num_T):
            args['T'] = T_lst[i]
            Cor_lst[i] = np.array([0 for _ in range(num_k)], dtype=float)
            for j in range(n):
                if method == 'MH':
                    result = Metropolis(args)
                elif method == 'Wolff':
                    result = Wolff(args)
                Cor_lst[i] += result[3]
            Cor_lst[i] = np.divide(Cor_lst[i], n)

        print(Cor_lst)

        coef = []
        for i in range(len(T_lst)):
            T = T_lst[i]
            cor_lst = np.array(Cor_lst[i])
            cor_lst = -np.log(cor_lst)
            regr = linear_model.LinearRegression()
            regr.fit(k_lst.reshape(-1, 1), cor_lst)
            print('Score: %s' % regr.score(k_lst.reshape(-1, 1), cor_lst))
            print('Coefficient of %s : %s' % (round(T,2),regr.coef_))
            coef.append(1/regr.coef_)
            if T > critical:
                plt.plot(k_lst, cor_lst, label = 'T = %s' % round(T,2), linewidth = 0, marker = 'o', markersize = 8)
                plt.plot(k_lst, regr.predict(k_lst.reshape(-1, 1)), color='red', linewidth=2)
        plt.xlabel(r'$k$')
        plt.ylabel(r'$-\log(\Gamma(k))$')
        plt.legend()
        plt.show()

        print(coef)

        coef = np.array(coef)
        plt.plot(T_lst, coef, linewidth = 1, marker = 'o', markersize = 8)
        plt.xlabel(r'$T$')
        plt.ylabel(r'$\xi$')
        plt.show()

        for i in range(num_T):
            if T_lst[i] <= critical:
                i += 1
            else:
                break

        T_lst = T_lst[i+1:]
        coef = coef[i+1:]
        log_coef = np.log(np.abs(coef))
        log_epsilon = np.log(np.abs(T_lst / critical - np.ones(len(T_lst))))


        regr = linear_model.LinearRegression()
        regr.fit(log_epsilon.reshape(-1, 1), log_coef)
        print('Score: %s' % regr.score(log_epsilon.reshape(-1, 1), log_coef))
        print('Coefficient: %s' % regr.coef_)

        plt.subplot(1,2,1)
        plt.plot(log_epsilon, log_coef, linewidth = 0, marker = 'o', markersize = 8)
        plt.plot(log_epsilon, regr.predict(log_epsilon.reshape(-1, 1)), color='red', linewidth=2)
        plt.xlabel(r'$\log \epsilon$')
        plt.ylabel(r'$\log \xi$')
        plt.title('Regression')


        plt.subplot(1,2,2)
        plt.plot(T_lst, coef, linewidth = 0, marker = 'o', markersize = 8)
        plt.plot(T_lst, np.exp(regr.predict(log_epsilon.reshape(-1, 1))), color='red', linewidth=2)
        plt.xlabel('$T$')
        plt.ylabel(r'$\xi$')
        plt.title('Prediction')
        plt.show()
    elif choice[1:] == 'belowCor':
        assert method == 'Wolff'
        if choice[0] == '2':
            N = 24; num_k = 10; num_T = 21; T_1 = 0.9; T_2 = 0.98
        else:
            N = 16; num_k = 7; num_T = 11; T_1 = 1.75; T_2 = 1.8
        assert T_2 < critical
        args = parameter_specify(q, N, d)
        k_lst = np.array([i for i in range(1, num_k+1)])
        args['J'] = J; args['kB'] = kB; args['n'] = n; args['k_lst'] = k_lst; args['h'] = 0; args['critical'] = critical
        T_lst = np.linspace(T_1, T_2, num_T)
        Cor_lst = [[0 for _ in range(num_k)] for _ in  range(num_T)]
        for i in range(num_T):
            args['T'] = T_lst[i]
            result = Wolff(args)
            Cor_lst[i] = result[3]

        print(Cor_lst)

        def curve(x, a, xi, c):
            return a * np.exp(-x/xi) + c

        coef = []
        for i in range(len(T_lst)):
            T = T_lst[i]
            cor_lst = np.array(Cor_lst[i])
            popt, pcov = optimize.curve_fit(curve, k_lst, cor_lst)
            a, xi, c = popt
            coef.append(xi)
            std_a, std_xi, std_c = np.sqrt(np.diag(pcov))
            print('T:%s' % T)
            print('Regression a:%s, xi:%s, c:%s' %(a,xi,c))
            print('Std: a:%s, xi:%s, c:%s'%(std_a, std_xi, std_c))
            enhanced_k_lst = np.linspace(1, num_k, 100)
            predicted_cor_lst = curve(enhanced_k_lst, *popt)
            plt.plot(k_lst, cor_lst, label = 'T = %s' % round(T,3), linewidth = 0, marker = 'o', markersize = 8)
            plt.plot(enhanced_k_lst, predicted_cor_lst, color='red', linewidth=2)
        plt.xlabel(r'k')
        plt.ylabel(r'$\Gamma(k)$')
        plt.legend()
        plt.show()

        log_coef = np.log(np.abs(coef))
        log_epsilon = np.log(np.abs(T_lst / critical - np.ones(len(T_lst))))

        regr = linear_model.LinearRegression()
        regr.fit(log_epsilon.reshape(-1, 1), log_coef)
        print('Score: %s' % regr.score(log_epsilon.reshape(-1, 1), log_coef))
        print('Coefficient: %s' % regr.coef_)

        plt.subplot(1,2,1)
        plt.plot(log_epsilon, log_coef, linewidth = 0, marker = 'o', markersize = 8)
        plt.plot(log_epsilon, regr.predict(log_epsilon.reshape(-1, 1)), color='red', linewidth=2)
        plt.xlabel(r'$\log \epsilon$')
        plt.ylabel(r'$\log \xi$')
        plt.title('Regression')


        plt.subplot(1,2,2)
        plt.plot(T_lst, coef, linewidth = 0, marker = 'o', markersize = 8)
        plt.plot(T_lst, np.exp(regr.predict(log_epsilon.reshape(-1, 1))), color='red', linewidth=2)
        plt.xlabel('$T$')
        plt.ylabel(r'$\xi$')
        plt.title('Prediction')
        plt.show()










if __name__ == '__main__':
    # main_simulation('2Mag', 'MH')
    # main_simulation('2Mag', 'MH')
    # main_simulation('2Hal', 'MH')
    # main_simulation('3Hal', 'MH')
    # main_simulation('2Cor', 'MH')
    # main_simulation('3Cor', 'MH')
    # main_simulation('2Mag', 'Wolff')
    # main_simulation('2Hal', 'Wolff')
    # main_simulation('2Mag','Wolff')
    # main_simulation('3Hal', 'Wolff')
    # main_simulation('2Cor', 'Wolff')
    # main_simulation('3Cor', 'Wolff')
    # main_simulation('2belowCor', 'Wolff')
    # main_simulation('3belowCor', 'Wolff')