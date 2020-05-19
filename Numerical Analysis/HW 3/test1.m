clear all
syms w1 w2 w3 x2 x3 y1
f = [w1+w2+w3-1/2;
    w2*x2+w3*x3-1/6;
    w1*y1+w3*(1-x3)-1/6;
    w2*x2^2+w3*x3^2-1/12;
    w1*y1^2+w3*(1-x3)^2-1/12;
    w3*x3*(1-x3)-1/24];
df = jacobian(f,[w1 w2 w3 x2 x3 y1]);
f0 = matlabFunction(f);
df0 = matlabFunction(df);
f1 = @(x) f0(x(1),x(2),x(3),x(4),x(5),x(6));
df1 = @(x) df0(x(1),x(2),x(3),x(4),x(5),x(6));

opts = []; opts.flag = 0;
epsilon = 0.3;
x0 = [1/6+rand()*epsilon,1/6+rand()*epsilon,1/6+rand()*epsilon,1/2+rand()*epsilon,1/2+rand()*epsilon,1/2+rand()*epsilon]';
x1 = newton_homotopy(f1, df1, x0, opts);
f1(x1)
