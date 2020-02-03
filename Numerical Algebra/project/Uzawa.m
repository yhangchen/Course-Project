function [U,V,P] = Uzawa(U, V, P, F, G, A, B1, B2, D, N, v, alpha)
% implemention of exact Uzawa iteration v times, by built-in linear solver.
    for i = 1:v
		F0 = F - B1 * P; G0 = G - B2 * P;
		U = A\F0; V = A\G0;
		P = P + alpha * ((B1'*U + B2'*V)-D);
    end
end