function x = least_square(A, b)
A0 = A'*A; b0 = A'*b;
x = QR_lin(A0,b0);
end
