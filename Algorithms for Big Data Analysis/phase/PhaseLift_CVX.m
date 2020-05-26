function [x,err] = PhaseLift_CVX(A, At, y, x0, lambda, epsilon)

% min_X Tr(X) + lambda*||X||_1
%  s.t. ||sum_i b(i) - Tr(A(i,:)'*A(i,:)*X)||< epsilon,
%       X >= 0.        
%
% and returns the eigenvector with the largest 
% eigenvalue of X. Handles both complex and real A and 
% Xs. 

if nargin < 6
    epsilon = 1e-3;
end
if nargin < 5
    lambda = 0;
end

%  We need direct access to the entires of A, 
%  so convert a function handle to a dense matrix

N = length(y);
n = length(x0);

if ~isnumeric(A)
    A0 = zeros(N,n);
    for i = 1:n
        ei = zeros(n,1); ei(i) = 1;
        A0(:,i) = A(ei);
    end
end
A = A0;



% [N, n]=size(A);

% call CVX to solve the convex problem



cvx_begin sdp quiet
    cvx_solver mosek
    variable X(n,n) hermitian
    minimize norm(diag(A*X*A')-y,1) + lambda*norm( X(:), 1 )
    subject to
        X >= 0;
cvx_end


% pick out the vector associated with the largest singular value
if any(isnan(X))
    x=NaN(size(X,1),1);
    disp('NaN returned')
else
    [u,s,v]=svds(X,1);
    x=u*sqrt(s);
end

err = norm(x0 - exp(-1i*angle(x0'*x)) * x)/norm(x0);

end
