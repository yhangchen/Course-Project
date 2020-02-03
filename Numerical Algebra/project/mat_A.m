function A = mat_A(N)
R = mat_R(N-1);
T = mat_T(N-1);
A = kron(R,speye(N-1)) + kron(speye(N),T);
end
