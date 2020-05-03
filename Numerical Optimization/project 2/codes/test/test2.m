function [f, g] = test2(x, n)
%{

    f(x)=\sum_{j=1}^{n/4}[
                (x_{4j-3}+10x_{4j-2})^2+5(x_{4j-1}-x_{4j})^2+(x_{4j-2}-2x_{4j-1})^4+10(x_{4j-3}-x_{4j})^4
            ]

 Parameters:

   Input, integer N, the number of variables.

   Input, real X(N), the argument of the objective function.

   Output, real F, the value of the objective function.

   Output, real G(N), the gradient of the objective function.

from:
https://people.sc.fsu.edu/~jburkardt/m_src/test_opt/test_opt.html
%}

if mod(n, 4); error('n should be a multiple of 4'); end
  f = 0.0;

  for j = 1 : 4 : n

    if ( j + 1 <= n )
      xjp1 = x(j+1);
    else
      xjp1 = 0.0;
    end

    if ( j + 2 <= n )
      xjp2 = x(j+2);
    else
      xjp2 = 0.0;
    end

    if ( j + 3 <= n )
      xjp3 = x(j+3);
    else
      xjp3 = 0.0;
    end

    f1 = x(j) + 10.0 * xjp1;

    if ( j + 1 <= n )
      f2 = xjp2 - xjp3;
    else
      f2 = 0.0;
    end

    if ( j + 2 <= n )
      f3 = xjp1 - 2.0 * xjp2;
    else
      f3 = 0.0;
    end

    if ( j + 3 <= n )
      f4 = x(j) - xjp3;
    else
      f4 = 0.0;
    end

    f = f +        f1 * f1 ...
          +  5.0 * f2 * f2 ...
          +        f3 * f3 * f3 * f3 ...
          + 10.0 * f4 * f4 * f4 * f4;

  end
  g = zeros ( n, 1 );

  for j = 1 : 4 : n

    if ( j + 1 <= n )
      xjp1 = x(j+1);
    else
      xjp1 = 0.0;
    end

    if ( j + 2 <= n )
      xjp2 = x(j+2);
    else
      xjp2 = 0.0;
    end

    if ( j + 3 <= n )
      xjp3 = x(j+3);
    else
      xjp3 = 0.0;
    end

    f1 = x(j) + 10.0 * xjp1;
    df1dxj = 1.0;
    df1dxjp1 = 10.0;

    if ( j + 1 <= n )
      f2 = xjp2 - xjp3;
      df2dxjp2 = 1.0;
      df2dxjp3 = -1.0;
    else
      f2 = 0.0;
      df2dxjp2 = 0.0;
      df2dxjp3 = 0.0;
    end

    if ( j + 2 <= n )
      f3 = xjp1 - 2.0 * xjp2;
      df3dxjp1 = 1.0;
      df3dxjp2 = -2.0;
    else
      f3 = 0.0;
      df3dxjp1 = 0.0;
      df3dxjp2 = 0.0;
    end

    if ( j + 3 <= n )
      f4 = x(j) - xjp3;
      df4dxj = 1.0;
      df4dxjp3 = -1.0;
    else
      f4 = 0.0;
      df4dxj = 0.0;
      df4dxjp3 = 0.0;
    end

    g(j) = 2.0 * f1 * df1dxj + 10.0 * 4.0 * f4^3 * df4dxj;

    if ( j+1 <= n )
      g(j+1) = 2.0 * f1 * df1dxjp1 + 4.0 * f3^3 * df3dxjp1;
    end

    if ( j+2 <= n )
      g(j+2) = 2.0 * 5.0 * f2 * df2dxjp2 + 4.0 * f3^3 * df3dxjp2;
    end

    if ( j+3 <= n )
      g(j+3) = 2.0 * 5.0 * f2 * df2dxjp3 + 10.0 * 4.0 * f4^3 * df4dxjp3;
    end

  end
end
