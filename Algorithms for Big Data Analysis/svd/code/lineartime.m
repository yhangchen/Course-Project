function [err, sval, svec] = lineartime(A, k, c, p)
% Input: 
% A: the matrix
% k: the rank of the approximation
% p: the prob vector
% 1<= k <=c <=n
%
% Return:
% err: ||A-QQ'A||
% sing_val: the k largest approximate singular values of A
% sing_vec: k singular vectors corresponding to the singular values in
% 'sing_val'

[m, n] = size(A); err = 0;
sval = 0; svec = 0;
C = zeros(m,c);
ind =  randsample(1:n,c,true,p);
parfor i = 1:c
    C(:,i) = A(:,ind(i))/sqrt(c*p(ind(i)));
end
[V,Sigma,~] = svd(C'*C);
Sigma = sqrt(diag(Sigma));
V = V(:,1:k);
sval = sval + Sigma(1:k);
svec = svec + C*V*diag(1./sval);
A_ap = svec'*A; A_ap = svec*A_ap;
err = err + norm(A-A_ap);

% sval = sval/10;
% svec = svec/10;
% err = err/10;
end
