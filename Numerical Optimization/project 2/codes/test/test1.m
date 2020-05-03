function [f, g] = test1(x, n)
%{

    f(x)=\sum_{i=1}^n \left\{
        n+i-\sum_{j=1}^n \left[
            \delta_{ij}\sin(x_j)+(i\delta_{ij}+1)\cos(x_j)
        \right]
        \right\}^2,\quad n=1,2,\cdots


 Parameters:

   Input, integer N, the number of variables.

   Input, real X(N), the argument of the objective function.

   Output, real F, the value of the objective function.

   Output, real G(N), the gradient of the objective function.

from:
https://people.sc.fsu.edu/~jburkardt/m_src/test_opt/test_opt.html
%}
  s1 = sum ( cos ( x(1:n) ) );

  f = 0.0;
  for j = 1 : n
    t =  ( n + j ) - sin ( x(j) ) - s1 - j * cos ( x(j) );
    f = f + t * t;
  end

  g = zeros ( n, 1 );


  s2 = 0.0;
  for j = 1 : n
    t =  ( n + j ) - sin ( x(j) ) - s1 - j * cos ( x(j) );
    s2 = s2 + t;
    g(j) = ( j * sin ( x(j) ) - cos ( x(j) ) ) * t;
  end

  for j = 1 : n
    g(j) = 2.0 * ( g(j) + sin ( x(j) ) * s2 );
  end

end


