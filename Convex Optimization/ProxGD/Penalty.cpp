#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

#include <cfloat> // use DBL_EPSILON to avoid zero division error.
#include <stdarg.h>
#include "ProxGD.h"
#include <string>
using namespace Eigen;
using namespace std;

Penalty::Penalty(string mode0, int n, ...)
{
    /* The first input is the mode of the penalty function.

    For general functions, the second input is penalty parameter mu, the third input is double alpha (for elastic net), double matrix D (for General Group Lasso). If we use the TV norm penalty, we need to provide the nrows and ncols of the input x to generate sparse matrix D_sp, and call GLasso solver. We store D_T, D_sp_T, being the transpose of D and D_sp for future use.

    For indicator functions, the second input is int (only for L_0 and rank constraints). For box constriants, we need to specify the upper and lower bound matrix, U and L separately. For half space, we need to provide the matrix and constant; for affine functions, we vectorize matrix x, and provide matrix A (A.ncols=x.size()) and vector b.  
    
    */
    mode = mode0;
    if (mode.substr(0, 3) == "Ind")
    {
        if (n >= 1)
        {
            va_list args;
            va_start(args, n);
            if ((mode == "Ind_rank") || (mode == "Ind_L_0"))
            {
                R0 = va_arg(args, int);
            }
            else if (mode == "Ind_box")
            {
                L = va_arg(args, MatrixXd);
                U = va_arg(args, MatrixXd);
            }
            else if (mode == "Ind_affine")
            {
                A = va_arg(args, MatrixXd);
                b = va_arg(args, MatrixXd);
                AAT = A * A.transpose(); // store AAT in advance for future use.
            }
            else if (mode == "Ind_half")
            {
                A = va_arg(args, MatrixXd);
                constant = va_arg(args, double);
            }
            else
            {
                R = va_arg(args, double);
            }
        }
    }
    else
    {
        va_list args;
        va_start(args, n);

        mu = va_arg(args, double);
        if (mode == "Elastic")
            alpha = va_arg(args, double);
        else if (mode == "GLasso")
        {
            D_T = va_arg(args, MatrixXd);
            D_T = D_T.transpose();
        }
        else if (mode == "TV_1D")
        {
            int nrows = va_arg(args, int); // # of rows of input vector.
            typedef Eigen::Triplet<int> T;
            std::vector<T> tripletList;
            tripletList.reserve((2 * nrows - 1) * sizeof(int));
            for (int i = 0; i < nrows - 1; i++)
            {
                tripletList.push_back(T(i, i, 1));
                tripletList.push_back(T(i, i + 1, -1));
            }
            Eigen::SparseMatrix<double> mat(nrows - 1, nrows);
            mat.setFromTriplets(tripletList.begin(), tripletList.end());
            D_sp_T = mat.transpose();
        }
        else if (mode == "TV_2D")
        {
            int nrows = va_arg(args, int); // # of rows of input vector.
            int ncols = va_arg(args, int); // # of cols of input vector.
            typedef Eigen::Triplet<int> T;
            std::vector<T> tripletList;
            tripletList.reserve(2 * (2 * ncols * nrows - nrows - ncols) * sizeof(int));
            for (int i = 0; i < nrows - 1; i++)
            {
                for (int j = 0; j < ncols; j++)
                {
                    tripletList.push_back(T(j * (nrows - 1) + i, j * (nrows) + i, 1));
                    tripletList.push_back(T(j * (nrows - 1) + i, j * (nrows) + i + 1, -1));
                }
            }
            int tmp = ncols * (nrows - 1);
            for (int i = 0; i < nrows * (ncols - 1); i++)
            {
                tripletList.push_back(T(tmp + i, i, 1));
                tripletList.push_back(T(tmp + i, i + nrows, 1));
            }
            Eigen::SparseMatrix<double> mat(2 * ncols * nrows - nrows - ncols, ncols * nrows);
            mat.setFromTriplets(tripletList.begin(), tripletList.end());
            D_sp_T = mat.transpose();
        }
        va_end(args);
    }
}

double Penalty::get_mu()
{
    return mu;
}
string Penalty::get_mode()
{
    return mode;
}

double Penalty::h(MatrixXd x) // wrapper
{
    if (mode == "L_0")
    {
        return L_0(x);
    }
    else if (mode == "L_1")
        return L_1(x);
    else if (mode == "L_2")
        return L_2(x);
    else if (mode == "L_12")
        return L_12(x);
    else if (mode == "L_21")
        return L_21(x);
    else if (mode == "L_inf")
        return L_inf(x);
    else if (mode == "Nuclear")
        return Nuclear(x);
    else if (mode == "Log_barrier")
        return Log_barrier(x);
    else if (mode == "Elastic")
        return Elastic(x);
    else if (mode == "GLasso")
        return GLasso(x);
    else if (mode == "TV_1D")
        return TV_1D(x);
    else if (mode == "TV_2D")
        return TV_2D(x);

    else
        throw "incorrect objective function.";
}

MatrixXd Penalty::prox_h(MatrixXd x, double new_mu) // wrapper.
{
    mu = new_mu;
    if (mode == "L_0")
        return L_0_prox(x);
    else if (mode == "L_1")
        return L_1_prox(x);
    else if (mode == "L_2")
        return L_2_prox(x);
    else if (mode == "L_12")
        return L_12_prox(x);
    else if (mode == "L_21")
        return L_21_prox(x);
    else if (mode == "L_inf")
        return L_inf_prox(x);
    else if (mode == "Nuclear")
        return Nuclear_prox(x);
    else if (mode == "Log_barrier")
        return Log_barrier_prox(x);
    else if (mode == "Elastic")
        return Elastic_prox(x);
    else if (mode == "GLasso")
        return GLasso_prox(x);
    else if (mode == "TV_1D")
        return TV_1D_prox(x);
    else if (mode == "TV_2D")
        return TV_2D_prox(x);
    else if (mode == "Ind_L_0")
        return Ind_L_0_prox(x);
    else if (mode == "Ind_L_1")
        return Ind_L_1_prox(x);
    else if (mode == "Ind_L_F")
        return Ind_L_F_prox(x);
    else if (mode == "Ind_L_inf")
        return Ind_L_inf_prox(x);
    else if (mode == "Ind_L_inf_2")
        return Ind_L_inf_2_prox(x);
    else if (mode == "Ind_box")
        return Ind_box_prox(x);
    else if (mode == "Ind_positive")
        return Ind_positive_prox(x);
    else if (mode == "Ind_negative")
        return Ind_negative_prox(x);
    else if (mode == "Ind_half")
        return Ind_half_prox(x);
    else if (mode == "Ind_affine")
        return Ind_affine_prox(x);
    else if (mode == "Ind_nuclear")
        return Ind_nuclear_prox(x);
    else if (mode == "Ind_psd")
        return Ind_psd_prox(x);
    else if (mode == "Ind_rank")
        return Ind_rank_prox(x);
    else
        throw "incorrect Penalty function.";
}

MatrixXd Penalty::L1_soft(MatrixXd x, double mu0)
{
    return (x.array().sign() * (x.array().abs() - mu0).max(0)).matrix();
}

MatrixXd Penalty::Ind_L_0_prox(MatrixXd x)
{
    MatrixXd x1 = x.cwiseAbs();
    Map<VectorXd> v(x1.data(), x1.size());
    std::sort(v.data(), v.data() + v.size());
    double threshold = v(v.size() - R0) - DBL_EPSILON;
    ArrayXXd mask = (x.cwiseAbs().array() > threshold).cast<double>();
    return x.cwiseProduct(mask.matrix());
}

MatrixXd Penalty::Ind_L_1_prox(MatrixXd x)
{
    double l1 = x.cwiseAbs().sum();
    double linf = x.cwiseAbs().maxCoeff();
    if (l1 < R)
    {
        return x;
    }
    double u = linf;
    double l = 0;
    double m = (u + l) / 2.0;
    double vm = (Penalty::L1_soft(x, m)).cwiseAbs().sum();
    while (fabs(vm - R) > FLT_EPSILON) //bisection method.
    {
        if (vm > R)
        {
            l = m;
        }
        else
        {
            u = m;
        }
        m = (u + l) / 2.0;
        vm = (Penalty::L1_soft(x, m)).cwiseAbs().sum();
    }
    return Penalty::L1_soft(x, m);
}

MatrixXd Penalty::Ind_L_F_prox(MatrixXd x)
{
    return x * min(1.0, R / (x.norm() + DBL_EPSILON));
}

MatrixXd Penalty::Ind_L_inf_prox(MatrixXd x)
{
    return x.cwiseMax(-R).cwiseMin(R);
}
MatrixXd Penalty::Ind_L_inf_2_prox(MatrixXd x)
{
    int n = x.rows();
    for (int i = 0; i < n; i++)
        x.row(i) *= min(1.0, R / (x.row(i).norm() + DBL_EPSILON));
    return x;
}

MatrixXd Penalty::Ind_box_prox(MatrixXd x)
{
    assert(L.rows() == x.rows());
    assert(L.cols() == x.cols());
    assert(U.rows() == x.rows());
    assert(U.cols() == x.cols());
    return x.cwiseMax(L).cwiseMin(U);
}

MatrixXd Penalty::Ind_positive_prox(MatrixXd x)
{
    return x.cwiseMax(0);
}

MatrixXd Penalty::Ind_negative_prox(MatrixXd x)
{
    return x.cwiseMin(0);
}

MatrixXd Penalty::Ind_half_prox(MatrixXd x)
{
    assert(A.rows() == x.rows());
    assert(A.cols() == x.cols());
    Map<VectorXd> x1(x.data(), x.size());
    Map<VectorXd> a(A.data(), A.size());
    return x - A / (A.norm() + DBL_EPSILON) * (x1.dot(a) - constant);
}

MatrixXd Penalty::Ind_affine_prox(MatrixXd x)
{
    int n = x.rows(), r = x.cols(), m = A.rows();
    assert(A.cols() == n * r);
    assert(b.rows() == A.rows());
    Map<VectorXd> v(x.data(), x.size());
    VectorXd result = A * v - b;
    result = AAT.ldlt().solve(result);
    result = v - A.transpose() * result;
    Map<MatrixXd> mat_result(result.data(), x.rows(), x.cols());
    return mat_result;
}

MatrixXd Penalty::Ind_nuclear_prox(MatrixXd x)
{
    MatrixXd result = x;
    JacobiSVD<MatrixXd> svd(x, ComputeThinU | ComputeThinV);
    result = svd.matrixU() * Penalty::Ind_L_1_prox(svd.singularValues()).asDiagonal() * svd.matrixV().transpose();
    return result;
}

MatrixXd Penalty::Ind_psd_prox(MatrixXd x)
{
    assert(x.cols() == x.rows());
    MatrixXd result = x;
    JacobiSVD<MatrixXd> svd((x + x.transpose()) / 2.0, ComputeThinU | ComputeThinV);
    result = svd.matrixU() * svd.singularValues().cwiseMax(0).asDiagonal() * svd.matrixV().transpose();
    return result;
}

MatrixXd Penalty::Ind_rank_prox(MatrixXd x)
{
    JacobiSVD<MatrixXd> svd((x + x.transpose()) / 2.0, ComputeThinU | ComputeThinV);
    return svd.matrixU().leftCols(R0) * svd.singularValues().head(R0).asDiagonal() * svd.matrixV().leftCols(R0).transpose();
}

double Penalty::L_12(MatrixXd x)
{
    double result = 0.0;
    for (int i = 0; i < x.rows(); i++)
    {
        result += x.row(i).norm();
    }
    return result;
}

double Penalty::L_21(MatrixXd x)
{
    double result = 0.0;
    for (int j = 0; j < x.cols(); j++)
    {
        result += x.col(j).norm();
    }
    return result;
}

double Penalty::Nuclear(MatrixXd x)
{
    JacobiSVD<MatrixXd> svd(x, ComputeThinU | ComputeThinV);
    double result = svd.singularValues().sum();
    return result;
}

double Penalty::L_inf(MatrixXd x)
{
    assert(x.cols() == 1);
    return x.lpNorm<Infinity>();
}

double Penalty::L_0(MatrixXd x)
{
    double result = 0.0;
    assert(x.cols() == 1);
    for (int i = 0; i < x.rows(); i++)
    {
        result += (fabs(x(i, 0)) > DBL_EPSILON) ? 1 : 0;
    }
    return result;
}

double Penalty::GLasso(MatrixXd x)
{
    MatrixXd y;
    if (D_T.size() > 0) // dense
    {
        assert(D_T.rows() == x.rows());
        y = D_T.transpose() * x;
    }
    else // sparse
    {
        assert(D_sp_T.rows() == x.rows());
        y = D_sp_T.transpose() * x;
    }

    double result = 0.0;
    for (int i = 0; i < y.rows(); i++)
    {
        result += y.row(i).norm();
    }
    return result;
}

double Penalty::Log_barrier(MatrixXd x)
{
    return -x.array().log().sum();
}

int Penalty::is_positive(MatrixXd x)
{
    for (int i = 0; i < x.rows(); i++)
    {
        for (int j = 0; j < x.cols(); j++)
        {
            if (x(i, j) <= DBL_EPSILON)
            {
                throw "x must be positive.";
            }
        }
    }
    return 0;
}

double Penalty::TV_1D(MatrixXd x)
{
    assert(x.cols() == 1);
    return Penalty::GLasso(x);
}

double Penalty::TV_2D(MatrixXd x)
{
    Map<VectorXd> v(x.data(), x.size());
    return Penalty::GLasso(v);
}

MatrixXd Penalty::TV_1D_prox(MatrixXd x)
{
    assert(x.cols() == 1);
    return Penalty::GLasso_prox(x);
}

MatrixXd Penalty::TV_2D_prox(MatrixXd x)
{
    Map<VectorXd> v(x.data(), x.size());
    MatrixXd result = Penalty::GLasso_prox(v);
    Map<MatrixXd> result_new(result.data(), x.rows(), x.cols());
    return result_new;
}

MatrixXd Penalty::L_12_prox(MatrixXd x)
{
    MatrixXd result = x;
    for (int i = 0; i < x.rows(); i++)
    {
        if (x.row(i).norm() < mu)
        {
            result.row(i) *= 0;
        }
        double flag = 1 - mu / x.row(i).norm();
        result.row(i) *= max(flag, 0.0);
    }
    return result;
}

MatrixXd Penalty::L_21_prox(MatrixXd x)
{
    MatrixXd result = x;
    for (int i = 0; i < x.cols(); i++)
    {
        if (x.col(i).norm() < DBL_EPSILON)
        {
            result.col(i) *= 0;
        }
        double flag = 1 - mu / x.col(i).norm();
        result.col(i) *= max(flag, 0.0);
    }
    return result;
}

MatrixXd Penalty::Nuclear_prox(MatrixXd x)
{
    MatrixXd result = x;
    JacobiSVD<MatrixXd> svd(x, ComputeThinU | ComputeThinV);
    result = svd.matrixU() * (svd.singularValues().array() - mu).max(0).matrix().asDiagonal() * svd.matrixV().transpose();
    return result;
}

MatrixXd Penalty::L_0_prox(MatrixXd x)
{
    MatrixXd result = x;
    assert(x.cols() == 1);
    double sqrt_mu = sqrt(2 * mu);
    for (int i = 0; i < x.rows(); i++)
    {
        result(i, 0) = (fabs(x(i, 0)) > sqrt_mu) ? x(i, 0) : 0;
    }
    return result;
}

MatrixXd Penalty::L_inf_prox(MatrixXd x)
{
    assert(x.cols() == 1);
    MatrixXd result = x;
    MatrixXd x_abs(x.size() + 1, 1);
    x_abs << x.cwiseAbs(), 0;
    sort(x_abs.data(), x_abs.data() + x_abs.size(), greater<double>());
    // rearrange from highest to lowest, add zeros in the end
    double t = 0;
    // see the report for determining t.
    double cum = 0;
    int i = 0;
    for (i = 0; (i < x_abs.size() - 1) && (cum < mu); i++)
    {
        cum += (i + 1) * (x_abs(i) - x_abs(i + 1));
    }
    // t in [x_abs(i - 1), x_abs(i)]
    if (cum < mu)
    {
        return 0 * result; //
    }
    t = x_abs(i) + (cum - mu) / i;
    // obtain t
    for (int i = 0; i < x.size(); i++)
    {
        if (x(i) > t)
        {
            result(i, 0) = t;
        }
        else if (x(i) < -t)
        {
            result(i, 0) = -t;
        }
    }
    return result;
}

MatrixXd Penalty::GLasso_prox(MatrixXd x)
{
    MatrixXd mx = -x;
    MatrixXd *x_p = &mx;
    if (D_T.size() > 0) // dense.
    {
        assert(D_T.rows() == x.rows());
        MatrixXd W_init = MatrixXd::Zero(D_T.cols(), x.cols());
        Objective f_obj("Frob", &D_T, x_p);
        Penalty h_penalty("Ind_L_inf_2", 1, mu);
        Result W = ProxGD_one_step(f_obj, h_penalty, "BB", &D_T, x_p, W_init, mu);
        MatrixXd W0 = W.min_point();
        return x + D_T * W0;
    }
    else // sparse.
    {
        assert(D_sp_T.rows() == x.rows());
        MatrixXd W_init = MatrixXd::Zero(D_sp_T.cols(), x.cols());
        Objective_Sparse f_obj("Frob", &D_sp_T, x_p);
        Penalty h_penalty("Ind_L_inf_2", 1, mu);
        Result W = ProxGD_Sparse_one_step(f_obj, h_penalty, "BB", &D_sp_T, x_p, W_init, mu);
        MatrixXd W0 = W.min_point();
        return x + D_sp_T * W0;
    }
}

MatrixXd Penalty::Log_barrier_prox(MatrixXd x)
{
    ArrayXXd xray = x.array();
    return (xray + (xray * xray + 4 * mu).sqrt()).matrix() / 2.0;
}

double Penalty::L_1(MatrixXd x)
{
    assert(x.cols() == 1);
    return Penalty::L_12(x);
}

double Penalty::L_2(MatrixXd x)
{
    assert(x.cols() == 1);
    return Penalty::L_21(x);
}

double Penalty::Elastic(MatrixXd x)
{
    assert((0.0 <= alpha) && (alpha <= 1.0));
    return Penalty::L_12(x) * alpha + (1 - alpha) * x.squaredNorm() / 2.0;
}

MatrixXd Penalty::Elastic_prox(MatrixXd x)
{
    assert((0.0 <= alpha) && (alpha <= 1.0));
    for (int i = 0; i < x.rows(); i++)
    {
        double normx = x.row(i).norm();
        if (normx <= alpha * mu)
        {
            x.row(i) *= 0.0;
        }
        double flag = (1 - alpha * mu / normx) / (1.0 + (1.0 - alpha) * mu);
        x.row(i) *= max(flag, 0.0);
    }
    return x;
}

MatrixXd Penalty::L_1_prox(MatrixXd x)
{
    assert(x.cols() == 1);
    return Penalty::L_12_prox(x);
}

MatrixXd Penalty::L_2_prox(MatrixXd x)
{
    assert(x.cols() == 1);
    return Penalty::L_21_prox(x);
}