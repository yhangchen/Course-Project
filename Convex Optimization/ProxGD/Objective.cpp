#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "ProxGD.h"
#include <math.h>
#include <assert.h>
using namespace Eigen;
using namespace std;

Objective::Objective(string mode, MatrixXd *A, MatrixXd *b) : mode(mode), A(A), b(b){};
// initialize the class, including mode ("Frob" or "Log"), matrix A, b.

MatrixXd *Objective::get_A()
{
	return A;
}

MatrixXd *Objective::get_b()
{
	return b;
}
string Objective::get_mode()
{
	return mode;
}

double Objective::f(MatrixXd x)
{
	if (mode == "Frob")
		return Frob(x);
	else if (mode == "Log")
		return Logistic(x);
	else
		throw "Unknown mode!";
}
MatrixXd Objective::grad_f(MatrixXd x)
{
	if (mode == "Frob")
		return Frob_grad(x);
	else if (mode == "Log")
		return Logistic_grad(x);
	else
		throw "Unknown mode!";
};

double Objective::Frob(MatrixXd x)
{
	return ((*A) * x - (*b)).squaredNorm() / 2.0;
}

double Objective::Logistic(MatrixXd x)
{
	assert(x.cols() == 1);
	return (1 + (-((*A) * x).cwiseProduct(*b).array()).exp()).log().sum() / (double)x.rows();
}

MatrixXd Objective::Frob_grad(MatrixXd x)
{
	return (*A).transpose() * (*A * x - *b);
}

MatrixXd Objective::Logistic_grad(MatrixXd x)
{
	int m = x.rows();
	MatrixXd result = x, Ax = *A * x;
	ArrayXd tmp = (*b).array();
	tmp = (-Ax.array() * tmp).exp();
	tmp *= (*b).array() / (1.0 + tmp) / (double)m;
	return -(*A).transpose() * tmp.matrix();
}

int Objective::check(MatrixXd x) // check the input data.
{
	int m = (*A).rows(),
		r = (*A).cols(), r1 = x.rows(), n = x.cols(), m1 = (*b).rows(), n1 = (*b).cols();
	if ((m != m1) || (r != r1) || (n != n1))
		throw "Incompatible size.";
	return 0;
}