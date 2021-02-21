#include <iostream>
#include <stdarg.h>
#include <stdlib.h>
#include <random>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "ProxGD.h"
#include <assert.h>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>
#include <cfloat>
#include <string>
using namespace Eigen;
using namespace std;

Result::Result(int iter, MatrixXd x, double min_value, double penalty_value) : x(x), iter(iter), min_value(min_value), penalty_value(penalty_value), t(-1), exact_x_in(false){};
// This is the constructor. It can initialize every memeber variable.

void Result::show()
{
	cout << "iterations: " << iter << endl;
	cout << "min value: " << min_value << endl;
	cout << "penalty value: " << penalty_value << endl;
	if (t >= 0)
	{
		cout << "CPU time: " << t << endl;
	}
	if (exact_x_in)
	{
		cout << "error: " << err_function(x, exact_x) << endl;
	}
	cout << "sparsity: " << sparsity(x) << endl;
	// This function can show each member variable.
}

MatrixXd Result::min_point()
{
	return x; // This function can return the optimal point.
}

double Result::min_loss()
{
	return min_value; // This function can return the optimal value.
}

int Result::iterations()
{
	return iter; // This function can return the time of iterations.
}

int Result::modify_iter(int k)
{
	iter = k;
	return k; // Modify the iter to k.
}

double Result::add_time(double time)
{
	// set the cpu time;
	t = time;
	return t;
}

void Result::add_exact_x(MatrixXd x)
{
	// Add the exact solution.
	exact_x_in = 1;
	exact_x = x;
}

double err_function(MatrixXd x, MatrixXd x0)
{
	// Calculate the error.
	MatrixXd delta_x = x - x0;
	return delta_x.norm() / (1 + max(x0.norm(), x.norm()));
}

double sparsity(MatrixXd x)
{
	// To calculate the sparsity.
	int m = x.rows();
	int n = x.cols();
	int sp = 0;
	int threshold = x.cwiseAbs().maxCoeff() * 1e-6;
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			sp += (fabs(x(i, j)) > threshold) ? 1 : 0;
		}
	}
	return sp * 1.0 / (m * n);
}
