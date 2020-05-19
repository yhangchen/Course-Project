clear all
syms x1 y1 w2
f = [(x1^2+(1-x1-y1)^2+y1^2)*(1/6-w2/3)+w2/9-1/12;
    (x1^2*y1+(1-x1-y1)^2*x1+x1^2*(1-x1-y1))*(1/6-w2/3)+w2/27-1/60;
    (x1^3+(1-x1-y1)^3+y1^3)*(1/6-w2/3)+w2/27-1/20];
df = jacobian(f,[w2 x1 y1]);

f0 = matlabFunction(f);
df0 = matlabFunction(df);

f1 = @(x) f0(x(1),x(2),x(3));
df1 = @(x) df0(x(1),x(2),x(3));

opts = []; opts.flag = 0;
epsilon = 0.5;
x0 = zeros(3,1);
x1 = newton_homotopy(f1, df1, x0, opts);
x1
f1(x1)

