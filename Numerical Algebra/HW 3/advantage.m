k = 1e7;%¿Éµ÷Õû
A = [k/sqrt(2)+1/2 -k/sqrt(6)+sqrt(3)/2
    -k/sqrt(2) k/sqrt(6)
    k/sqrt(2)-1/2 -k/sqrt(6)-sqrt(3)/2];
b = [3 2 -1]';
x = least_square(A,b);
%QR
% x = cho(A'*A,A'*b);
%NC
% x = gauss_max(A'*A,A'*b);
%NG
r = norm(A*x-b)
err = norm(x-[1,sqrt(3)]')
