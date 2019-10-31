function b = improved_cho(A,b)%improved_cholesky decomposition
    n = length(A);
    v = zeros(n,1);
    for j = 1:n
        for i = 1:j-1
            v(i)=A(j,i)*A(i,i);
        end
        A(j,j)=A(j,j)-A(j,1:j-1)*v(1:j-1);
        A(j+1:n,j)=(A(j+1:n,j)-A(j+1:n,1:j-1)*v(1:j-1))/A(j,j);
    end
    for j = 1:n-1
        b(j+1:n)=b(j+1:n)-b(j)*A(j+1:n,j);
    end
    for i = 1:n
        b(i)=b(i)/A(i,i);
    end
    for j = n:-1:2
        b(1:j-1)=b(1:j-1)-b(j)*A(j,1:j-1)';
    end
end
