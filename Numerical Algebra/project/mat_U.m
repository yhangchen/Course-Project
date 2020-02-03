function out = mat_U(N)
u = zeros(N-1,N);
for i = 1:N-1
    for j = 1:N
        u(i,j) = (1-cos(2*pi*i/N))*sin(2*pi*(j-0.5)/N);
    end
end
v = -u;
U = u(:); V = v(:);
out = [U;V];
end