#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "ProxGD.h"
#include <math.h>
#include <assert.h>
using namespace Eigen;
using namespace std;

Objective_Sparse::Objective_Sparse(string mode, SparseMatrix<double> *A, MatrixXd *b) : mode(mode), A(A), b(b){};

SparseMatrix<double> *Objective_Sparse::get_A()
{
    return A;
}

MatrixXd *Objective_Sparse::get_b()
{
    return b;
}
string Objective_Sparse::get_mode()
{
    return mode;
}

double Objective_Sparse::f(MatrixXd x)
{
    if (mode == "Frob")
        return Frob(x);
    else if (mode == "Log")
        return Logistic(x);
    else
        throw "Unknown mode!";
}
MatrixXd Objective_Sparse::grad_f(MatrixXd x)
{
    if (mode == "Frob")
        return Frob_grad(x);
    else if (mode == "Log")
        return Logistic_grad(x);
    else
        throw "Unknown mode!";
};

double Objective_Sparse::Frob(MatrixXd x)
{
    return ((*A) * x - (*b)).squaredNorm() / 2.0;
}

double Objective_Sparse::Logistic(MatrixXd x)
{
    assert(x.cols() == 1);
    return (1 + (-(*A * x).cwiseProduct(*b).array()).exp()).log().sum() / (double)x.rows();
}

MatrixXd Objective_Sparse::Frob_grad(MatrixXd x)
{
    return (*A).transpose() * (*A * x - *b);
}

MatrixXd Objective_Sparse::Logistic_grad(MatrixXd x)
{
    int m = x.rows();
    MatrixXd result = x, Ax = *A * x;
    ArrayXd tmp = (*b).array();
    tmp = (-Ax.array() * tmp).exp();
    tmp *= (*b).array() / (1.0 + tmp) / (double)m;
    return -(*A).transpose() * tmp.matrix();
}

int Objective_Sparse::check(MatrixXd x)
{
    int m = (*A).rows(),
        r = (*A).cols(), r1 = x.rows(), n = x.cols(), m1 = (*b).rows(), n1 = (*b).cols();
    if ((m != m1) || (r != r1) || (n != n1))
        throw "Incompatible size.";
    return 0;
}
