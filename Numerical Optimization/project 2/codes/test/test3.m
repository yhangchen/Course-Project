function [f, g] = test3(x, n)
%{

    f(x)=\sum_{i=2}^n[i(2x_i-x_{i-1})^2]

 Parameters:

   Input, integer N, the number of variables.

   Input, real X(N), the argument of the objective function.

   Output, real F, the value of the objective function.

   Output, real G(N), the gradient of the objective function.
%}

sub_d = [2:n]*(-2);
d = 5*[2:(n-1)]+1;
d = [2 d 4*n];
A = sparse([1:n-1 2:n 1:n], [2:n, 1:n-1 1:n], [sub_d sub_d d]);
f = x' * ( A * x );
g = 2 * A * x;
end

