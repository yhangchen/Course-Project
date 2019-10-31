function b = cho(A,b)%cholesky decomposition
    n = length(A);
    for k = 1:n
        A(k,k)=sqrt(A(k,k));
        A(k+1:n,k)=A(k+1:n,k)/A(k,k);
        for j = k+1:n
            A(j:n,j)=A(j:n,j)-A(j:n,k)*A(j,k);
        end
    end
    for j = 1:n-1
        b(j)=b(j)/A(j,j);
        b(j+1:n)=b(j+1:n)-b(j)*A(j+1:n,j);
    end
    for j = n:-1:2
        b(j)=b(j)/A(j,j);
        b(1:j-1)=b(1:j-1)-b(j)*A(j,1:j-1)';
    end
end