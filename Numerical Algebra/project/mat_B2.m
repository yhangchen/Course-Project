function B2 = mat_B2(N)
x = zeros(N*(N-1),1); y1 = zeros(N*(N-1),1); 
y2 = zeros(N*(N-1),1); one = ones(N*(N-1),1);
for i = 1:N-1
    for j = 1:N
        x((i-1)*N+j) = (j-1)*(N-1)+i;
        y1((i-1)*N+j) = j + (i-1)*N;
        y2((i-1)*N+j) = j + i*N;
    end
end
x = [x;x]; y = [y1;y2];
B2 = sparse(x,y,[-one;one],N*(N-1),N^2,2*N*(N-1));
end

        