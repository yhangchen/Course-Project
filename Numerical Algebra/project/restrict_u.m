function R = restrict_u(N)
% restriction matrix of velocity.
R = 1/8*kron(mat_R1(N/2),mat_R2(N/2));
end