#include <iostream>
#include <stdarg.h>
#include <Eigen/Dense>
#include "ProxGD.h"
#include <assert.h>
#include <algorithm>
using namespace Eigen;
using namespace std;

double line_search(Objective &f_obj, string tmode, MatrixXd x, double gamma, int n, ...)
{
	// This is the function of line searching, to decide the step size t. f_obj is the goal function.
	// tmode is the mode of line searching. x is the starting point. It will return t.
	// gamma is the parameter of the Armijo rule. n is the number of the variable parameters.
	// The variable parameters are parameters of Armijo and BB.

	double t;

	if (tmode == "Armijo")
	{
		// Armijo rule

		va_list args;
		va_start(args, n);
		double fh = va_arg(args, double);
		va_end(args);

		double l_t = 0;
		double r_t = 1;
		t = 1;
		bool flag = 1; // When flag = 1, we will extend the section of line searching.
		bool cond;	   // cond denotes whether the Armijo rule.

		while (r_t - l_t > 1e-8)
		{
			// If the section of search is very small, then stop.

			cond = (f_obj.f(x - t * f_obj.grad_f(x)) < f_obj.f(x) - gamma * t * f_obj.grad_f(x).squaredNorm()); // Judge whether the Armijo rule is satisfied.

			if (cond)
			{
				// If it is satisfied, decide whether to extend the section of line searching.
				if (flag)
				{
					l_t = r_t;
					r_t *= 2;
					t = r_t;
					// If needed, extend the section of line searching.
				}
				else
				{
					break;
					// Otherwise break.
				}
			}
			else
			{
				// If it is not satisfied, start narrowing the section of line searching.

				flag = 0; // Set the flag as 0 to stopping extending the section.

				r_t = t;
				t = (r_t + l_t) / 2; // Narrow the section and modify the step size.
			}
		}
	}
	else if (tmode == "BB")
	{
		// BB

		va_list args;
		va_start(args, n);
		double *fhs = va_arg(args, double *);
		int M = va_arg(args, int);
		MatrixXd delta_x = va_arg(args, MatrixXd);
		va_end(args);
		MatrixXd delta_g = f_obj.grad_f(x) - f_obj.grad_f(x - delta_x);
		t = delta_x.squaredNorm() / (delta_x.array().cwiseProduct(delta_g.array())).sum();
		while (f_obj.f(x - t * f_obj.grad_f(x)) > *(max_element(fhs, fhs + M)) - gamma * t * f_obj.grad_f(x).squaredNorm())
		{
			// If it is not satisfied, reduce the step size.

			t /= 2;
		}
	}
	return t;
}