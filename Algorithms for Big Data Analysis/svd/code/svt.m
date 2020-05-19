function out = svt(R, M, mu, tol, range)
% matrix size m,n
% R: observed location 0-1 matrix, size m*n
% M: observed data sparse matrix,
% k_max: maximum of iter
% mu: parameter
% epsilon: stopping criterion
% l: increment.
% range: elementwise range of M, usually used in graphic reconstruction.
[n1, n2] = size(R);
lb = range(1);
ub = range(2);
M = R.*M;
m = sum(sum(R));
delta = 1.2*n1*n2/m;
if mu == 0
    mu = 5 * n2;
end
l = 5;
normFM = norm(R.*M,'fro');
norm2M = svds(M, 1);
k0 = ceil(mu/delta/norm2M);
Y = k0*delta*M;
Error = []; Norm = []; out = []; Rank = [];
r = 0;
for k = 1:1000
    s = r + 1;
    delta = delta/(1+5e-4);
    [U, D, V] = svds(Y, s);
    while D(s, s) > mu
        s = s + l;
        [U, D, V] = svds(Y, s);
    end
    for i = s:-1:1
        if D(i,i) > mu
            break
        end
    end
    r = max(i, r);
    d = diag(D(1:r,1:r));
    d = d - mu;
    U = U(:,1:r); V = V(:,1:r);
    X = min(max((U*diag(d))*V',lb),ub);
    diff = R.*(M-X);
    stopping = norm(diff,'fro')/normFM;
    Error = [Error stopping];
    Rank = [Rank r];
    Norm = [Norm sum(abs(d))];
    fprintf('iter: '+string(k)+', rank: '+string(r)+' error: '+ string(stopping)+'\n');
    if stopping <= tol
        break
    end
    Y = Y+delta*diff;
end
out.X = X; out.Norm = Norm; out.Train = Error; out.Rank = Rank; 
end