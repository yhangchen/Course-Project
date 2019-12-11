function err = err_estimate(m)
%input size of grid m(even division)
%output error(L2 norm)
X = 0:1/m:1;
mu = finite_element(@(x) 1, X);
err = 1/12;
for i = 1:m
    err = err + (mu(i+1)-mu(i))^2*m - (mu(i+1)-mu(i))*(1-(2*i-1)/m);
end
err = sqrt(err);
end
