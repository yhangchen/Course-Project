clear all
syms x y z w2

f = [(1/24-w2/4)*(1-x-y)*(x+y)+w2/16-1/120;
    (1/24-w2/4)*(1-x-z)*(x+z)+w2/16-1/120;
    (1/24-w2/4)*(1-z-y)*(z+y)+w2/16-1/120;
    (1/24-w2/4)*(x^3+y^3+z^3+(1-x-y-z)^3)+w2/64-1/120];
    
df = jacobian(f, [w2, x, y, z]);
f0 = matlabFunction(f);
df0 = matlabFunction(df);
f1 = @(x) f0(x(1),x(2),x(3),x(4));
df1 = @(x) df0(x(1),x(2),x(3),x(4));

opts = []; opts.flag = 0;
x0 = zeros(4,1);
[x1,out] = newton_homotopy(f1, df1, x0, opts)
