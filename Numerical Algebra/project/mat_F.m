function F = mat_F(N)
F0 = zeros(N-1,N);
for i = 1:N-1
    for j = 1:N
        x = i/N; y = (j-0.5)/N;
        F0(i,j) = fun_f(x,y);
    end
end
e1 = zeros(N,1); e1(1) = 1;
eN = zeros(N,1); eN(N) = 1;
base = zeros(N-1,1);
for i = 1:N-1
    base(i) = i/N;
end
vec_b = fun_b(base); vec_t = fun_t(base);
F = F0(:) + N*(kron(e1,vec_b)+kron(eN,vec_t));
end

function out = fun_f(x,y)
out = -4*pi^2*(2*cos(2*pi*x)-1)*sin(2*pi*y)+x^2;
end

function out = fun_b(x)
out = -2*pi*(1-cos(2*pi*x));
end

function out = fun_t(x)
out = 2*pi*(1-cos(2*pi*x));
end