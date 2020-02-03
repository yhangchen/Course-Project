function G = mat_G(N)
G0 = zeros(N,N-1);
for i = 1:N
    for j = 1:N-1
        x = (i-0.5)/N; y = j/N;
        G0(i,j) = fun_g(x,y);
    end
end
e1 = zeros(N,1); e1(1) = 1;
eN = zeros(N,1); eN(N) = 1;
base = zeros(N-1,1);
for i = 1:N-1
    base(i) = i/N;
end
vec_l = fun_l(base); vec_r = fun_r(base);
G0t = G0';
G = G0t(:) + N*(kron(e1,vec_l)+kron(eN,vec_r));
end

function out = fun_g(x,y)
out = 4*pi^2*(2*cos(2*pi*y)-1)*sin(2*pi*x);
end

function out = fun_l(y)
out = 2*pi*(1-cos(2*pi*y));
end

function out = fun_r(y)
out = -2*pi*(1-cos(2*pi*y));
end