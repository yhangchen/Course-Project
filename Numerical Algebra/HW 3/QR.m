function [A,d] = QR(A)
[m,n]=size(A);d = zeros(n,1);
for j = 1:n
    if j<m
        [v,beta]=house(A(j:m,j));
        w = v'*A(j:m,j:n);
        A(j:m,j:n)=A(j:m,j:n)-v*beta*w;
        d(j)=beta;
        A(j+1:m,j)=v(2:m-j+1);
    end
end
end

