A = [10 1 2 3 4; 1 9 -1 2 -3; 2 -1 7 3 -5; 3 2 3 12 -1; 4 -3 -5 -1 15];
b = [12; -27; 14; -17; 12];
real = [1, -2, 3, -2, 1]';
A_eigenvalue = eig(A)
[x1, k1, t1] = Gauss_Seidel(A, b, 1e-15, inf);
[x2, k2, t2] = Jacobi(A, b, 1e-15, inf);
[x3, k3, t3] = conjugate_gradient(A, b, 1e-6, inf);
e1 = norm(x1 - real); e2 = norm(x2 - real); e3 = norm(x3 - real);
e1
e2
e3
k1
k2
k3
t1
t2
t3