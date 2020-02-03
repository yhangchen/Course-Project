function BtB = mat_BtB(N)
R = mat_R(N-1);
BtB = kron(speye(N),R)+kron(R,speye(N));
end