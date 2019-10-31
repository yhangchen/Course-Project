function b = gauss(A,b,bool)
%解方程Ax=b，bool='T'表示使用列主元Gauss消去，bool='F'表示不选主元。
n = length(A);
x0 = A\b;%真解。
if strcmp(bool,'F')
    for k = 1:n-1
        if A(k,k) ~=0
            A(k+1:n,k)=A(k+1:n,k)/A(k,k);
            A(k+1:n,k+1:n)=A(k+1:n,k+1:n)-A(k+1:n,k)*A(k,k+1:n);
        end
        if A(k,k)==0
            ME = MException('Singular');
            throw(ME);
        end
    end
    %解下三角
    for j = 1:n-1
        b(j+1:n)=b(j+1:n)-b(j)*A(j+1:n,j);
    end
    %解上三角
    for j = n:-1:2
        b(j)=b(j)/A(j,j);
        b(1:j-1)=b(1:j-1)-b(j)*A(1:j-1,j);
    end
    b(1)=b(1)/A(1,1);
    err = norm(b-x0)
end
u = zeros(n-1,1);
if strcmp(bool,'T')
    for k = 1:n-1
        maximum = abs(A(k,k));
        max_ind = k;
        for k1 = k+1:n
            if abs(A(k1,k))> maximum
                maximum = abs(A(k1,k));
                max_ind = k1;
            end
        end
        u(k) = max_ind;
        A([k,u(k)],:)=A([u(k),k],:);
        b([k,u(k)])=b([u(k),k]);
        if A(k,k) ~=0
            A(k+1:n,k)=A(k+1:n,k)/A(k,k);
            A(k+1:n,k+1:n)=A(k+1:n,k+1:n)-A(k+1:n,k)*A(k,k+1:n);
        end
        if A(k,k)==0
            ME = MException('Singular');
            throw(ME);
        end
    end
    for j = 1:n-1
        b(j+1:n)=b(j+1:n)-b(j)*A(j+1:n,j);
    end
    for j = n:-1:2
        b(j)=b(j)/A(j,j);
        b(1:j-1)=b(1:j-1)-b(j)*A(1:j-1,j);
    end
    b(1)=b(1)/A(1,1);
    err = norm(b-x0)
end
end
    