function [f, g, G] = test_disc_ie(x)
%{
    implement the discrete integral equation
        minimize f(x)=\sum_{i=1}^m f_i^2(x), x\in \mathbb{R}^n
        where m = n;
        f_i(x)&
    =x_i+h[(1-t_i)\sum_{j=1}^i t_j(x_j+t_j+1)^3+t_i\sum_{j=i+1}^n (1-t_j)(x_j+t_j+1)^3]/2

    input:
        x\in \mathbb{R}^n,
    output:
       f: value; g: gradient; G: Hessian
%}
n = length(x); % x must be size nx1
h = 1/(n+1);
x1 = x + [1:n]'*h + 1;
x2 = x1.^2;
x3 = x2.*x1;
t = [1:n]'*h;
sum1 = 0;
sum2 = sum((1-t) .* x3);
f_r = zeros(n, 1); J = zeros(n, n); S = zeros(n, n);
for i = 1:n
    ti = t(i);
    sum1 = sum1 + ti * x3(i);
    sum2 = sum2 - (1-ti) * x3(i);
    fi = x(i) + 0.5 * h * ((1-ti) * sum1 + ti * sum2);    
    f_r(i) = fi;
    
    for j = 1:n
        if j < i
            J(i,j) = 3 * h * (1-ti) * t(j) * (x(j) + t(j) + 1)^2 / 2;
        elseif j == i
            J(i,j) = 1 + 3 * h * (1-ti) * t(j) * (x(j) + t(j) + 1)^2 / 2;
        else
            J(i,j) = 3 * h * ti * (1-t(j)) * (x(j) + t(j) + 1)^2 / 2;
        end
    end
    
    dSi = zeros(n, 1);
    for j = 1:n
        if j <= i
            dSi(j) = 3 * h * (1-ti) * t(j) * (x(j) + t(j) + 1);
        else
            dSi(j) = 3 * h * ti * (1-t(j)) * (x(j) + t(j) + 1);
        end
    end
    
    S = S + fi * diag(dSi); % since \nabla^2 f_i(x) is diagnoal.    
end
f = sum(f_r.^2);
g = 2 * J' * f_r;
G = 2 * (J' * J + S);
end


