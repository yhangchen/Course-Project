#define _CRT_SECURE_NO_WARNINGS

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

void save_file(string name, string datatype, MatrixXd &A)
{
    // Save the matrix A as csv files. name is the file name without ".csv".
    // The file will be in the folder "data".
    stringstream ss;
    ss << "./results/" << name << datatype << ".csv ";
    ss >> name;
    ofstream fout;
    fout.open(name, ios::out | ios::trunc);
    for (int i = 0; i < A.rows(); i++)
    {
        for (int j = 0; j < A.cols(); j++)
        {
            fout << A(i, j) << ' ';
        }
        fout << '\n';
    }
    fout.close();
}

MatrixXd read_mat(int row, int col, string name)
{
    MatrixXd A(row, col);
    A.setZero();

    std::ifstream indata;
    indata.open(name);

    std::string line;
    std::stringstream lineStream;
    std::string cell;
    std::getline(indata, line);
    for (int i = 0; i < row; i++)
    {
        if (std::getline(indata, line))
        {
            lineStream.str("");
            lineStream.clear();
            lineStream << line;
            for (int j = 0; j < col; j++)
            {
                if (std::getline(lineStream, cell, ','))
                {
                    A(i, j) = std::stod(cell);
                    cell = "";
                }
                else
                {
                    cout << "This line is:" << endl;
                    cout << line << endl;
                    cout << "The matrix is not wide enough." << endl;
                    cout << "Now j = " << j << "." << endl;
                    cout << "Now i = " << i << "." << endl;
                    cout << "The name is: " << name << "." << endl;
                }
            }
        }
        else
        {
            cout << "The matrix is not high enough." << endl;
            cout << "Now i = " << i << "." << endl;
            cout << "The name is: " << name << "." << endl;
        }
    }

    indata.close();
    return A;
}

double denoise_mat(string name1, string name2, string fmode, string hmode, string tmode, int row, int col, double mu)
{
    string path1;
    string path2;
    string path3;
    SparseMatrix<double> A(row, row);
    MatrixXd min_point;
    A.setIdentity();
    cout << "Denoise Mat" << endl;
    for (int i = 0; i < 3; i++)
    {
        cout << "layer:" << i + 1 << endl;
        path1 = "./datasets/" + name1 + to_string(i + 1) + ".csv";
        MatrixXd x0 = read_mat(row, col, path1) / 255.0;
        path2 = "./datasets/" + name2 + to_string(i + 1) + ".csv";
        MatrixXd x_exact = read_mat(row, col, path2) / 255.0;
        Objective_Sparse f_obj(fmode, &A, &x0);
        Penalty h_penalty(hmode, 3, mu, row, col);
        Result res = ProxGD_Sparse_one_step(f_obj, h_penalty, tmode, &A, &x0, x0, mu);
        res.add_exact_x(x_exact);
        res.show();
        path3 = name2 + fmode + "_" + hmode + "_" + tmode + "_denoiseMat_" + to_string(i + 1) + ".csv";
        min_point = res.min_point();
        // save_file(path3, "_A", A);
        save_file(path3, "_b", x_exact);
        save_file(path3, "_x0", x0);
        save_file(path3, "_x", min_point);
    }

    return 0.0;
}

double denoise_vec(string name1, string name2, string fmode, string hmode, string tmode, int row, int col, double mu)
{
    string path1;
    string path2;
    string path3;
    stringstream ss;
    SparseMatrix<double> A(row, row);
    MatrixXd min_point;
    A.setIdentity();
    cout << "Denoise Vec" << endl;
    for (int i = 0; i < 2; i++)
    {
        if (i == 0)
            cout << "left channel" << endl;
        else
            cout << "right channel" << endl;
        path1 = "";
        path2 = "";
        path3 = "";
        ss.str("");
        ss.clear();
        ss << "./datasets/" << name1 << i + 1 << ".csv";
        ss >> path1;
        MatrixXd x0 = read_mat(row, 1, path1);
        ss.str("");
        ss.clear();
        ss << "./datasets/" << name2 << i + 1 << ".csv";
        ss >> path2;
        MatrixXd x_exact = read_mat(row, 1, path2);
        Objective_Sparse f_obj(fmode, &A, &x0);
        Penalty h_penalty(hmode, 3, mu, row, col);
        Result res = ProxGD_Sparse_one_step(f_obj, h_penalty, tmode, &A, &x0, x0, mu, 1e-6);
        res.add_exact_x(x_exact);
        res.show();
        path3 = name2 + fmode + "_" + hmode + "_" + tmode + "_denoiseVec_" + to_string(i + 1) + ".csv";
        min_point = res.min_point();
        save_file(path3, "_b", x_exact);
        save_file(path3, "_x0", x0);
        save_file(path3, "_x", min_point);
    }

    return 0.0;
}

void sparsity_test(Objective &f_obj, Penalty &h_penalty, string tmode, MatrixXd &exact_x)
{
    MatrixXd *A = f_obj.get_A();
    MatrixXd *b = f_obj.get_b();
    string fmode = f_obj.get_mode();
    string hmode = h_penalty.get_mode();
    double mu = h_penalty.get_mu();

    Eigen::MatrixXd x0((*A).cols(), (*b).cols());
    x0.setRandom();
    // Initialize A£¬b£¬x0, ecact_x.
    Result res = ProxGD(f_obj, h_penalty, tmode, A, b, x0, mu, 1e-10);
    res.add_exact_x(exact_x);
    res.show();
    // Calculate and show the result.
    string name = fmode + "_" + hmode + "_" + tmode;
    MatrixXd min_point = res.min_point();
    save_file(name, "_A", *A);
    save_file(name, "_b", *b);
    save_file(name, "_exact_x", exact_x);
    save_file(name, "_x0", x0);
    save_file(name, "_x", min_point);
}

void main_sparsity_vector_test()
{
    static default_random_engine e(0);
    int m = 1024;
    int n = 512;
    int r = 1;
    double mu = 1e-2;
    double sp = 1e-1; // sparsity.
    Eigen::MatrixXd A(m, n);
    A.setRandom();

    int *a = new int[n];
    for (int i = 0; i < n; i++)
    {
        a[i] = i;
    }
    std::shuffle(a, a + n, std::default_random_engine(0));
    // vector case
    Eigen::MatrixXd exact_x(n, r);
    exact_x.setZero();
    for (int i = 0; i < ceil(n * sp); i++)
    {
        exact_x.row(a[i]).setRandom();
    }
    Eigen::MatrixXd b = A * exact_x;

    Objective f_obj("Frob", &A, &b);
    string tmode = "BB";
    string hmode = "L_0";
    Penalty h_penalty(hmode, 1, mu);
    cout << hmode << endl;
    sparsity_test(f_obj, h_penalty, tmode, exact_x);
    hmode = "L_1";
    h_penalty = Penalty(hmode, 1, mu);
    cout << hmode << endl;
    sparsity_test(f_obj, h_penalty, tmode, exact_x);
    hmode = "L_2";
    h_penalty = Penalty(hmode, 1, mu);
    cout << hmode << endl;
    sparsity_test(f_obj, h_penalty, tmode, exact_x);
    hmode = "L_inf";
    h_penalty = Penalty(hmode, 1, mu);
    cout << hmode << endl;
    sparsity_test(f_obj, h_penalty, tmode, exact_x);
    hmode = "Elastic";
    double alpha = 0.5;
    h_penalty = Penalty(hmode, 2, mu, alpha);
    cout << hmode << endl;
    sparsity_test(f_obj, h_penalty, tmode, exact_x);
}

void main_sparsity_mat_test()
{
    // matrix case
    static default_random_engine e(0);
    int m = 1024;
    int n = 512;
    int r = 4;
    double mu = 1e-2;
    double sp = 1e-1; // sparsity.
    Eigen::MatrixXd A(m, n);
    A.setRandom();

    int *a = new int[n];
    for (int i = 0; i < n; i++)
    {
        a[i] = i;
    }
    std::shuffle(a, a + n, std::default_random_engine(0));
    // vector case
    Eigen::MatrixXd exact_x(n, r);
    exact_x.setZero();
    for (int i = 0; i < ceil(n * sp); i++)
    {
        exact_x.row(a[i]).setRandom();
    }
    Eigen::MatrixXd b = A * exact_x;

    Objective f_obj("Frob", &A, &b);
    string tmode = "BB";
    string hmode = "L_12";
    Penalty h_penalty(hmode, 1, mu);
    cout << hmode << endl;
    sparsity_test(f_obj, h_penalty, tmode, exact_x);
    hmode = "L_21";
    h_penalty = Penalty(hmode, 1, mu);
    cout << hmode << endl;
    sparsity_test(f_obj, h_penalty, tmode, exact_x);
    hmode = "Elastic";
    double alpha = 0.5;
    h_penalty = Penalty(hmode, 2, mu, alpha);
    cout << hmode << endl;
    sparsity_test(f_obj, h_penalty, tmode, exact_x);
}

int main()
{
    main_sparsity_vector_test(); // L0,L1,L2,Linf,Elastic
    main_sparsity_mat_test();    // L12,L21.
    denoise_vec("fragment2_noise", "fragment2", "Frob", "TV_1D", "BB", 52920, 1, 0.1);
    return 0;
}
