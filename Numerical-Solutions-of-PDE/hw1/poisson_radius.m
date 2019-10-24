function u = poisson_radius(f,g,m,n)
%m为r方向网格数，n为θ方向网格数。f，g定义与报告中相同。
A = sparse(m*(n-1),m*(n-1));
dr = 1/m; dth = pi/2/n;
%矩阵。
for k = 1:n-1
    for j = 2:m-1
        A((k-1)*m+j,(k-1)*m+j)=(j-0.5)*(2*j-1)*dth^2+2;
    end
end
for j = 2:m-1
    for k = 1:n-1
        A((k-1)*m+j,(k-1)*m+j+1) = -(j-0.5)*(j)*dth^2;
        A((k-1)*m+j,(k-1)*m+j-1) = -(j-0.5)*(j-1)*dth^2;
        if k > 1
            A((k-1)*m+j,(k-2)*m+j) = -1;
        end
        if k < n-1
            A((k-1)*m+j,(k)*m+j) = -1;
        end
    end
end

for k = 1:n-1
    A((k-1)*m+m,(k-1)*m+m)=(m-0.5)*dth^2*(3*m-1)+2;
    A((k-1)*m+m,(k-1)*m+m-1)=-(m-0.5)*dth^2*(m-1);
    A((k-1)*m+1,(k-1)*m+1)=2+dth^2/2;
    A((k-1)*m+1,(k-1)*m+2)=-dth^2/2;
    if k>1
       A((k-1)*m+m,(k-2)*m+m)=-1;
       A((k-1)*m+1,(k-2)*m+1)=-1;
    end
    if k<n-1
       A(k*m,(k+1)*m)=-1;
       A((k-1)*m+1,(k)*m+1)=-1;
    end
end
%右端项
b = zeros(m*(n-1),1);
for j = 1:m
    for k = 1:n-1
        b((k-1)*m+j)=f((j-0.5)*dr,k*dth)*dr^2*(j-0.5)^2*dth^2;
    end
end

for i = 1:m
    b(i)=b(i)+g((i-0.5)*dr,0);
    b(m*(n-2)+i)=b(m*(n-2)+i)+g((i-0.5)*dr,pi/2);
end
for i = 1:n-1
    b(m*i)=b(m*i)+m*(2*m-1)*dth^2*g(1,i*dth);
end
%u = A\b;
u=gauss_band(A,b,n);
end

function b = gauss_band(A,b,m)
%带宽2m+1
n = length(A);
for k = 1:n-1
    if A(k,k) ~=0
        s = min(k+m,n);
        A(k+1:s,k)=A(k+1:s,k)/A(k,k);
        A(k+1:s,k+1:s)=A(k+1:s,k+1:s)-A(k+1:s,k)*A(k,k+1:s);
    end
    if A(k,k)==0
        ME = MException('Singular');
        throw(ME);
    end
end
for j = 1:n-1
    s = min(j+m,n);
    b(j+1:s)=b(j+1:s)-b(j)*A(j+1:s,j);
end
for j = n:-1:2
    s = max(j-m,1);
    b(j)=b(j)/A(j,j);
    b(s:j-1)=b(s:j-1)-b(j)*A(s:j-1,j);
end
b(1)=b(1)/A(1,1);
end