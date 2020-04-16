function [W,flag] = model2(n)
% n size of the matrix
% condition number of output matrix, set to be n
% if flag ~= 1, the solution is not valid.
R = sprand(n,n,0.1);
[ii, jj]= find(R);
B = sparse(ii, jj, 0.5, n,n);
L = triu(B,1);
B = L + L';
f = @(x)cond(B+(x+5)*eye(n)) - n;
[delta,~,flag] = fsolve(f, 1);
W = B/delta + eye(n);
W = (W + W')/2;
try
    chol(W+eps*eye(n)); % check if positive definite
catch MSE
    flag = 0;
end
end
