function R = mat_R(N)
r1 = [1, 2 * ones(1, N-1), 1];
r2 = - ones(1, N);
R = sparse([1:N+1,1:N,2:N+1],[1:N+1,2:N+1,1:N], [r1, r2, r2]);
end