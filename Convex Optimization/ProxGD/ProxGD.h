#ifndef PROXGD
#define PROXGD
using namespace Eigen;
using namespace std;
#include <Eigen/Dense>
#include <Eigen/Sparse>

class Objective
{
public:
	Objective(string mode, MatrixXd *A, MatrixXd *b);
	double f(MatrixXd x);
	MatrixXd grad_f(MatrixXd x);
	double Frob(MatrixXd x);
	MatrixXd Frob_grad(MatrixXd x);
	double Logistic(MatrixXd x);
	MatrixXd Logistic_grad(MatrixXd x);
	int check(MatrixXd x);
	MatrixXd *get_A();
	MatrixXd *get_b();
	string get_mode();

private:
	string mode;
	MatrixXd *A;
	MatrixXd *b;
};

class Objective_Sparse
{
public:
	Objective_Sparse(string mode, SparseMatrix<double> *A, MatrixXd *b);
	double f(MatrixXd x);
	MatrixXd grad_f(MatrixXd x);
	double Frob(MatrixXd x);
	MatrixXd Frob_grad(MatrixXd x);
	double Logistic(MatrixXd x);
	MatrixXd Logistic_grad(MatrixXd x);
	int check(MatrixXd x);
	SparseMatrix<double> *get_A();
	MatrixXd *get_b();
	string get_mode();

private:
	string mode;
	SparseMatrix<double> *A;
	MatrixXd *b;
};

class Penalty
{
private:
	MatrixXd L1_soft(MatrixXd x, double mu0);
	double mu, alpha, R, constant;
	int R0;
	string mode;
	MatrixXd D_T, A, b, L, U, AAT;
	SparseMatrix<double> D_sp_T;

public:
	Penalty(string mode, int n, ...);
	double get_mu();
	string get_mode();
	double h(MatrixXd x);
	MatrixXd prox_h(MatrixXd x, double new_mu);
	int is_positive(MatrixXd x);
	double L_0(MatrixXd x);
	MatrixXd L_0_prox(MatrixXd x);
	double L_1(MatrixXd x);
	MatrixXd L_1_prox(MatrixXd x);
	double L_2(MatrixXd x);
	MatrixXd L_2_prox(MatrixXd x);
	double L_12(MatrixXd x);
	MatrixXd L_12_prox(MatrixXd x);
	double L_21(MatrixXd x);
	MatrixXd L_21_prox(MatrixXd x);
	double L_inf(MatrixXd x);
	MatrixXd L_inf_prox(MatrixXd x);
	double Nuclear(MatrixXd x);
	MatrixXd Nuclear_prox(MatrixXd x);
	double GLasso(MatrixXd x);
	MatrixXd GLasso_prox(MatrixXd x);
	double Log_barrier(MatrixXd x);
	MatrixXd Log_barrier_prox(MatrixXd x);
	double Elastic(MatrixXd x);
	MatrixXd Elastic_prox(MatrixXd x);
	double TV_1D(MatrixXd x);
	MatrixXd TV_1D_prox(MatrixXd x);
	double TV_2D(MatrixXd x);
	MatrixXd TV_2D_prox(MatrixXd x);

	MatrixXd Ind_L_0_prox(MatrixXd x);
	MatrixXd Ind_L_1_prox(MatrixXd x);
	MatrixXd Ind_L_F_prox(MatrixXd x);
	MatrixXd Ind_L_inf_prox(MatrixXd x);
	MatrixXd Ind_L_inf_2_prox(MatrixXd x);
	MatrixXd Ind_box_prox(MatrixXd x);
	MatrixXd Ind_positive_prox(MatrixXd x);
	MatrixXd Ind_negative_prox(MatrixXd x);
	MatrixXd Ind_half_prox(MatrixXd x);
	MatrixXd Ind_affine_prox(MatrixXd x);
	MatrixXd Ind_nuclear_prox(MatrixXd x);
	MatrixXd Ind_psd_prox(MatrixXd x);
	MatrixXd Ind_rank_prox(MatrixXd x);
};

class Result
{ // This class is used to store and show the result. See Result.cpp.
public:
	Result(int iter, MatrixXd x, double min_value, double penalty_value);
	void show(); // Show the result.
	MatrixXd min_point();
	double min_loss();
	int iterations();
	int modify_iter(int iter);
	double add_time(double t);
	void add_exact_x(MatrixXd x);

private:
	MatrixXd x;		  // Optimal x
	int iter;		  // Time of iterations
	double min_value; // Optimal value
	double t;		  // CPU time
	bool exact_x_in;  // When exact_x is added, this variable is true
	MatrixXd exact_x; // Exact solution
	double penalty_value;
};

Result ProxGD(Objective &f_obj, Penalty &h_penalty, string tmode, MatrixXd *A, MatrixXd *b, MatrixXd x0, double mu, double epsilon = 1e-6, double gamma = 0.5, int M = 2);
// To calculate the minimum of f(x)+g(x). See ProxGD.cpp.

Result ProxGD_one_step(Objective &f_obj, Penalty &h_penalty, string tmode, MatrixXd *A, MatrixXd *b, MatrixXd x0, double mu, double epsilon = 1e-6, double gamma = 0.5, int M = 2);
// To calculate the minimum of f(x)+g(x). See ProxGD.cpp.

Result ProxGD_Sparse(Objective_Sparse &f_obj, Penalty &h_penalty, string tmode, SparseMatrix<double> *A, MatrixXd *b, MatrixXd x0, double mu, double epsilon = 1e-6, double gamma = 0.5, int M = 2);
// To calculate the minimum of f(x)+g(x). See ProxGD.cpp.

Result ProxGD_Sparse_one_step(Objective_Sparse &f_obj, Penalty &h_penalty, string tmode, SparseMatrix<double> *A, MatrixXd *b, MatrixXd x0, double mu, double epsilon = 1e-6, double gamma = 0.5, int M = 2);
// To calculate the minimum of f(x)+g(x). See ProxGD.cpp.

double line_search(Objective &f_obj, string tmode, MatrixXd x, double gamma, int n, ...);
// The function of line searching. It can decide the step size t. See line_search.cpp.

double line_search_sparse(Objective_Sparse &f_obj, string tmode, MatrixXd x, double gamma, int n, ...);
// The function of line searching. It can decide the step size t. See line_search.cpp.

double err_function(MatrixXd x, MatrixXd x0);
// To calculate the error. See Result.cpp.

double sparsity(MatrixXd x);
// To calculate the sparsity. See Result.cpp.

#endif