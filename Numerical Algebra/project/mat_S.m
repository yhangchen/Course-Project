function S = mat_S(N)
s = ones(1, N);
S = -sparse([1:N, 1:N],[1:N,2:N+1],[s,-s]);
end