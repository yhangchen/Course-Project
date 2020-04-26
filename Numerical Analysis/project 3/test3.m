clear all
syms x y z
f = [(1-x-y)*(x+y)-1/5;
    (1-x-z)*(x+z)-1/5;
    (1-z-y)*(z+y)-1/5];
    
df = jacobian(f,[x,y,z]);
f0 = matlabFunction(f);
df0 = matlabFunction(df);
f1 = @(x) f0(x(1),x(2),x(3));
df1 = @(x) df0(x(1),x(2),x(3));


opts = []; opts.flag = 0;
x0 = rand(3,1);
x1 = newton_homotopy(f1, df1, x0, opts)
f1(x1)
