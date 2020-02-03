function R = restrict_p(N)
% restriction matrix of pressure.
R = 1/4*kron(mat_R1(N/2),mat_R1(N/2));
end
