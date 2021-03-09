function x = sor(A,b,w)
n = length(b);x = b;init_err=norm(A*b-b);err=init_err;
d = full(diag(A));l=full(diag(A,-1));u=full(diag(A,1));
%若用系数稀疏索引则速度较慢。
while err>init_err*10^(-10)
    x(1)=(1-w)*x(1)+w*(b(1)-u(1)*x(2))/d(1);
    for i = 2:n-1
        x(i)=(1-w)*x(i)+w*(b(i)-u(i)*x(i+1)-l(i-1)*x(i-1))/d(i);
    end
    x(n)=(1-w)*x(n)+w*(b(n)-l(n-1)*x(n-1))/d(n);
    err = norm(A*x-b);
end
end

