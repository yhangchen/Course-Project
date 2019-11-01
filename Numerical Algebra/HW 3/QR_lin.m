function x = QR_lin(A,b)
%解线性方程组
[m,n]=size(A);[A,d]=QR(A);
for j = 1:n
    if j<m
        beta = d(j);
        v = [1,A(j+1:m,j)']';
        w = v'*b(j:m);
        b(j:m)=b(j:m)-beta*v*w;
    end
end
for j = n:-1:2
    b(j)= b(j)/A(j,j);
    b(1:j-1)=b(1:j-1)-b(j)*A(1:j-1,j);
end
b(1)=b(1)/A(1,1);
x = b(1:n);
end

        