function [f, g, G] = test_brown_den(x,m)
%{
    implement the Brown_Dennis function
        minimize f(x)=\sum_{i=1}^m f_i^2(x), x\in \mathbb{R}^n
        where n=4, m\geq n
        f_i(x)=(x_1+t_i x_2-\exp(t_i))^2+(x_3+x_4\sin(t_i)-\cos(t_i))^2

    input:
        x\in \mathbb{R}^4, m
    output:
        out:
            out.f: value; out.g: gradient; out.G: Hessian
%}
n = 4;
if length(x) ~= n || m < n
    error('wrong input')
end

f = 0; g = zeros(n,1); G = zeros(n,n);

for i = 1:m
    ti = i / 5.0;
    f1 = x(1) + ti * x(2) - exp (ti);
    f2 = x(3) + sin (ti) * x(4) - cos (ti);
    
    f = f + (f1^2 + f2^2)^2;

    g(1) = g(1) + 4.0 * ( f1^3  + f1 * f2^2);
    g(2) = g(2) + 4.0 * ( f1^3 * ti + f1 * f2^2 * ti );
    g(3) = g(3) + 4.0 * ( f1^2 * f2 + f2^3 );
    g(4) = g(4) + 4.0 * ( f1^2 * f2 * sin(ti) + f2^3 * sin(ti));
    
    G(1,1) = G(1,1) + 12.0 * f1^2  + 4.0 * f2^2 ;
    G(1,2) = G(1,2) + 12.0 * f1^2 * ti + 4.0 * f2^2 * ti;
    G(1,3) = G(1,3) +  8.0 * f1 * f2 ;
    G(1,4) = G(1,4) +  8.0 * f1 * f2 * sin(ti);
    G(2,1) = G(2,1) + 12.0 * f1^2 * ti  +  4.0 * f2^2 * ti;
    G(2,2) = G(2,2) + 12.0 * f1^2 * ti * ti  +  4.0 * f2^2 * ti;
    G(2,3) = G(2,3) +  8.0 * f1 * f2 * ti;
    G(2,4) = G(2,4) +  8.0 * f1 * f2 * ti * sin(ti);
    G(3,1) = G(3,1) +  8.0 * f1 * f2;
    G(3,2) = G(3,2) +  8.0 * f1 * f2 * ti;
    G(3,3) = G(3,3) +  4.0 * f1^2 + 12.0 * f2^2;
    G(3,4) = G(3,4) +  4.0 * f1^2 * sin(ti) + 12.0 * f2^2 * sin(ti);
    G(4,1) = G(4,1) +  8.0 * f1 * f2 * sin(ti);
    G(4,2) = G(4,2) +  8.0 * f1 * f2 * sin(ti) * ti;
    G(4,3) = G(4,3) +  4.0 * f1^2 * sin(ti) + 12.0 * f2^2 * sin(ti);
    G(4,4) = G(4,4) +  4.0 * f1^2 * sin(ti) * sin(ti) + 12.0 * f2^2 * sin(ti) * sin(ti);
end
end