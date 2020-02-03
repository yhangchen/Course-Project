function R1 = mat_R1(N)
	I = speye(N);
	R1 = kron(I,[1,1]);
end