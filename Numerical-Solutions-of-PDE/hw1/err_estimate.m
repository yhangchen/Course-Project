function err = err_estimate(f,g,h,n)
%hÎªÕæ½â¡£
u = poisson(f,g,n);
v = zeros((n-1)^2,1);
for i = 1:n-1
    for j = 1:n-1
        v((i-1)*(n-1)+j)=h(i/n,j/n);
    end
end
err = norm(u-v,inf);
end