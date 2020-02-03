function [U, V, P] = inexact_Uzawa(U, V, P, F, G, A, B1, B2, D, N, v, alpha)
% implementation of inexact Uzawa as smoother.
    F0 = F - B1 * P; G0 = G - B2 * P;
    U = Gauss_Seidel(A, F0, U, v); V = Gauss_Seidel(A, G0, V, v);
    P = P + alpha * ((B1'*U + B2'*V)-D);
end