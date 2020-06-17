#include <cstdlib>
#include <iostream>
#include <windows.h>
#include <fstream>
#include <cstdio>
#include <random>
#include <cmath>
#include <list>
#include <algorithm>

using namespace std;

// generate double precision random number in [0,1]
inline double double_rand() {
	return (double(rand()) * (RAND_MAX + 1) + double(rand())) / (RAND_MAX + 1)/ (RAND_MAX + 1);
}

class Ising {
private:
	int N;
	int *state;
	double beta, J, kB;
    double _H = 0;

	// indexes of neighbours
    // right
	inline int r(int ind) {
		return (ind / N) * N + (ind % N + 1) % N;
	}
    // left
	inline int l(int ind) {
		return (ind / N) * N + (ind % N + N - 1) % N;
	}
    // down
	inline int d(int ind) {
		return (ind % N)  + ((ind / N + N - 1) % N) * N;
	}
    // upper
	inline int u(int ind) {
		return (ind % N)  + ((ind / N + 1) % N) * N;
	}


public:
    Ising(int _N, double _beta, double _J, double _kB) :
        N(_N), beta(_beta), J(_J), kB(_kB) {
    state = new int[N * N];
    init_state();
    Hamilton();
    }

    // random initiate the state.
    void init_state() {
        for (int i = 0; i < N * N; ++i)
            state[i] = 2 * (rand() % 2) - 1;
    }

    // calculate the Hamiltonian
    void Hamilton() {
        for (int i = 0; i < N * N; ++i){
            _H -= J / 2 * state[i] * (state[u(i)] + state[d(i)] + state[l(i)] + state[r(i)]);
        }
    }

	void flip(int i) {
        // update the Hamiltonian
        _H +=  + 2 * J * state[i] * (state[u(i)] + state[d(i)] + state[l(i)] + state[r(i)]);
        state[i] = -state[i];
	}

    // Metropolis algorithm
	void Metropolis() {
        int i = rand() % N;
        int j = rand() % N;
        int ind = i * N + j;
        double dH = 2 * J * state[ind] * (state[u(ind)] + state[d(ind)] + state[l(ind)] + state[r(ind)]);
        double p = min(1.0, exp(-beta * dH));
        double threshold = double_rand();
        if (threshold < p)
            flip(ind);
    }

	double H_out() {
		return _H;
	}
};


int main(){
	const int N = 100;
    const double kB = 1;
    const double J = 1;
    const double beta_opt = log(1.0 + sqrt(2.0)) / 2;
	const int cases = 31;
    const double d_cases = 31.0;
	double warming_up = 0;
    const double max_iter = 1e7;
	double betas[cases];
    double Us[cases];
    double Cs[cases];
    double H0;
    for (int count = 0; count < cases; ++count)
    {
        double beta = 0.1 + 0.05*count;
        betas[count] = beta;
        if (beta > beta_opt){
            	warming_up = 5*1e8;
        }
        else {
            	warming_up = 1e7;
        }
        // construct the model
        Ising model(N, beta, J, kB);
        // warmming up
        for (int j = 0; j < warming_up; ++j)
        {
            model.Metropolis();
        }
        double total_H = 0;
        double total_H2 = 0;
        // sampling
        for (int k = 0; k < max_iter; ++k)
        {
            model.Metropolis();
            H0 = model.H_out();
            total_H += H0;
            total_H2 += H0*H0;

        }
        double E_H = total_H / max_iter;
        double E_H2 = total_H2 / max_iter;
        // averaging
        Us[count] = E_H / N / N;
        Cs[count] = kB * beta * beta * (E_H2 - E_H * E_H) / N / N;
		cout << betas[count] << ", " << Us[count] << ", " << Cs[count]
				<< endl;
    }
	ofstream fout;
	fout.open("results.csv");
	for (int count = 0; count < cases; ++count)
		fout << betas[count] << ", " << Us[count] << ", " << Cs[count]
				<< endl;
	fout.close();
    system("PAUSE");
    return 0;
}
