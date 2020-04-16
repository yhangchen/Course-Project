function W = model1(n)
%n size of the matrix
W = zeros(n,n);
parfor i = 1:n
    for j = 1:n
        W(i,j)= 0.6^(abs(i-j));
    end
end
end