function [err, sval, svec] = prototype(A, r, q)
% Input: 
% A: the matrix
% r: the rank of the approximation
% q: the exponent
%
% Return:
% err: ||A-QQ'A||
% sing_val: the k largest approximate singular values of A
% sing_vec: k singular vectors corresponding to the singular values in
% 'sing_val'

[~, n] = size(A);
k = ceil(r/2);
% Stage A
Omega = randn(n,2*k);
Y = A * Omega; 
[Q , ~ ] = qr(Y);
for i = 1:q
Y = A' * Q;
[Q , ~ ] = qr(Y);
Y = A * Q;
[Q , ~ ] = qr(Y);
end

% Stage B
Q = Q(:,1:2*k);
B = Q' * A;  
[U0, Sigma, ~] = svd(B); 
U = Q * U0;

% output
err = norm(A-Q*B,'fro');
sval = diag(Sigma);
sval = sval(1:r);
svec = U(:,1:r);
end
