function out = mat_P(N)
p = zeros(N,N);
for i = 1:N
    for j = 1:N
        p(i,j) = ((i-0.5)/N)^3/3-1/12;
    end
end
out = p(:);
end