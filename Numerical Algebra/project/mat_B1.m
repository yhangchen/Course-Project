function B1 = mat_B1(N)
B1 = kron(speye(N), mat_S(N-1));
end