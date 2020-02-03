function [U, V, P] = DGS_1(U, V, P, F, G, A, B1, B2, BtB, D, N, v)
% DGS in Chen's paper
	D0 = spdiags(1./diag(BtB),0,N^2,N^2);
    G_S = tril(A);
    for j = 1:v
        U = U + G_S\(F - A*U - B1*P);
        V = V + G_S\(G - A*V - B2*P);
        dP = D0*(D-(B1'*U+B2'*V));
        U = U + B1*dP;
        V = V + B2*dP;
        P = P - BtB*dP;
    end
end