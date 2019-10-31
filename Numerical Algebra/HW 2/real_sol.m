function u = real_sol(N)
u = zeros((N-1)^2,1);
for i = 1:N-1
    for j = 1:N-1
        u((N-1)*(i-1)+j)=sin(2*pi*i/N)*sin(2*pi*j/N);
    end
end
end
