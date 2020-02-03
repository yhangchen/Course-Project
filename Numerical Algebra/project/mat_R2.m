function R2 = mat_R2(N)
	R1 = mat_R1(N-1);
	Z = sparse(N-1,1);
	R2 = [R1,Z]+[Z,R1];
end