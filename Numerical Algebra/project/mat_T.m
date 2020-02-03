function T = mat_T(N)
t1 = ones(1,N) * 2;
t2 = - ones(1,N-1);
T = sparse([1:N,1:N-1,2:N],[1:N,2:N,1:N-1],[t1, t2, t2]);
end

