#include <iostream>
#include <stdarg.h>
#include <stdlib.h>
#include <random>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "ProxGD.h"
#include <assert.h>
#include <string>
#include <ctime>
#include <cfloat>

using namespace Eigen;
using namespace std;

Result ProxGD(Objective &f_obj, Penalty &h_penalty, string tmode, MatrixXd *A, MatrixXd *b, MatrixXd x0, double mu, double epsilon, double gamma, int M)
{
	clock_t start, end;
	start = clock();
	// Start timing.

	int K = ceil(fabs(log2f(mu)));
	double muk = mu * pow(2, K);
	Result res(0, x0, 0, 0);
	int iter = 0;
	double epsilonk = epsilon * pow(2, K);
	for (int i = 0; i <= K; i++)
	{
		res = ProxGD_one_step(f_obj, h_penalty, tmode, A, b, res.min_point(), muk, epsilonk, gamma, M);
		muk /= 2;
		epsilonk /= 2;
		iter += res.iterations();
	}
	res.modify_iter(iter);
	end = clock();
	double endtime = (double)(end - start) / CLOCKS_PER_SEC;
	res.add_time(endtime); // Calculate the cpu time.
	return res;
}

Result ProxGD_one_step(Objective &f_obj, Penalty &h_penalty, string tmode, MatrixXd *A, MatrixXd *b, MatrixXd x0, double mu, double epsilon, double gamma, int M)
{
	// Calculate the min value of f(x)+g(x). Input fmode, hmode, tmode to decide f, h and the mode of line searching.
	// A and b are the paramters of f. x0 is the initial value of iteration. mu is the coefficent before h.
	// epsilon is the parameter of stopping rule. The returning value is a result class, containing the min value, the
	// optimal point and the times of iteration.

	clock_t start, end;
	start = clock();
	// Start timing.
	MatrixXd new_x = x0;
	MatrixXd prev_x = x0;
	MatrixXd delta_x = x0;
	MatrixXd x_star = x0;
	delta_x.setOnes();
	// Initialize variables. In each loop, both new_x and prev_x, denoting new x and old x respectively, will be used.
	// delta_x denotes the difference of new_x and prev_x. x_star is the middle value of the algorithm.
	// Initializing delta_x as an all-one matrix can guarantee the program to enter the while loop.

	int iter = 0;
	double t;
	double *fhs = new double[M];
	fhs[0] = f_obj.f(new_x);
	for (int i = 1; i < M; i++)
	{
		fhs[i] = fhs[0];
	}
	fhs[M - 1]++;
	// Initialize variables. t is the step size. iter is the time of iterations. fhs is the function value of the latter M
	// steps. Each step will update a function value. fhs[iter % M] is the current value and the latter one is the value of
	// last step and so on. The value before fhs[0] is fhs[M - 1]. Initialize fhs[0] and fhs[M - 1] to guarantee the program
	// to enter the while loop.

	while ((fabs(fhs[iter % M] - fhs[(iter + M - 1) % M]) > epsilon) && (delta_x.squaredNorm() > epsilon))
	{
		// when the change of x or f(x) + h(x) is less than epsilon, the loop will end.
		// Here we know the purpose of initializing fhs[M - 1] and delta_x.

		if (tmode == "Armijo")
		{
			t = line_search(f_obj, tmode, new_x, gamma, 1, (double)fhs[iter % M]);
		}
		else if (tmode == "BB")
		{
			t = line_search(f_obj, tmode, new_x, gamma, 3, fhs, M, delta_x);
		}
		// First, we need to get the step size by line searching

		x_star = new_x - t * f_obj.grad_f(new_x);
		// Then, we need to calculate x_star.

		prev_x = new_x;
		new_x = h_penalty.prox_h(x_star, t * mu); // Here, we use prox_h to calculate new x.

		delta_x = new_x - prev_x;
		// Update prev_x after updating new_x. Update delta_x at last.

		iter = iter + 1;
		// Update the iteration times.

		fhs[iter % M] = f_obj.f(new_x);
		// Update the function value.
	}
	double penalty_value;
	if (h_penalty.get_mode().substr(0, 3) != "Ind")
	{
		penalty_value = h_penalty.h(new_x);
	}
	else
	{
		penalty_value = 0.0;
	}

	Result res(iter, new_x, fhs[iter % M], penalty_value);

	end = clock();
	double endtime = (double)(end - start) / CLOCKS_PER_SEC;
	res.add_time(endtime); // Calculate the cpu time.
	return res;
}
