function [f, g] = test4(x, n, A)
%{
    \|B^2- A\|_2
 Parameters:

   Input, integer N, the number of variables.

   Input, real X(N), the argument of the objective function.

   Output, real F, the value of the objective function.

   Output, real G(N), the gradient of the objective function.
%}
if nargin < 3
    % set the default A
    xr = sin( [1:n^2].^2 );
    Br = reshape(xr, n, n)';
    A = Br * Br;
end
if length(x) ~= n^2 || size(A, 1) ~= n || size(A, 1) ~= n
    error('incorrect input'); 
end
B = reshape(x, n, n)';
f = sum(sum((B^2-A).^2));
C = B' * B';
g = 2 * (C * B + B * C - A' * B - B * A');
g = vec(g);
end
