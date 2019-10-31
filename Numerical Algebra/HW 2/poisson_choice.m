function [u,t] = poisson_choice(f,n,choice)
%构建迭代矩阵。
if strcmp(choice,'gauss')
    A = sparse(zeros((n-1)^2));
else
A = zeros((n-1)^2);
end

for i = 1:(n-1)^2
    A(i,i)=4;
end
for i = 1:(n-1)*(n-2)
    A(i,i+n-1)=-1;A(i+n-1,i)=-1;
end
for i = 1:(n-1)^2-1
    A(i,i+1)=-1;A(i+1,i)=-1;
end
for i = 1:n-2
    A((n-1)*i,(n-1)*i+1)=0;A((n-1)*i+1,(n-1)*i)=0;
end
b = zeros((n-1)^2,1);
for i = 1:n-1
    for j = 1:n-1
        b((i-1)*(n-1)+j)=f(i/n,j/n)/n^2;
    end
end
if strcmp(choice,'gauss')
    t0 = clock;
    u = gauss(A,b);
    t = etime(clock,t0);
end
if strcmp(choice,'gauss_band')
    t0 = clock;
    u = gauss_band(A,b,n);
    t = etime(clock,t0);
end
if strcmp(choice,'cholesky')
    t0 = clock;
    u = cho(A,b,n);
    t = etime(clock,t0);
end

%可视化解。
x = zeros((n-1)^2,1);
y = zeros((n-1)^2,1);
for i = 1:n-1
for j =1:n-1
x((i-1)*(n-1)+j)=i/n;
y((i-1)*(n-1)+j)=j/n;
end
end
% x_bd = zeros(4*n,1);y_bd = zeros(4*n,1);u_bd = zeros(4*n,1);
% for i = 1:n+1
%     x_bd(i)=0;
%     y_bd(i)=(i-1)/n;
%     x_bd(i+2*n)=1;
%     y_bd(i+2*n)=1-(i-1)/n;
% end
% for i =1:n-1
%     x_bd(i+n+1)=i/n;
%     y_bd(i+n+1)=1;
%     x_bd(i+3*n+1)=1-i/n;
%     y_bd(i+3*n+1)=0;
% end
% x = [x' x_bd']';y = [y' y_bd']';u = [u' u_bd']';
% scatter3(x,y,u);
end


function b = cho(A,b,m)%cholesky decomposition
    n = length(A);
    v = zeros(n,1);
    for j = 1:n
        s = max(j-m,1);
        u = min(j+m,n);
        for i = s:j-1
            v(i)=A(j,i)*A(i,i);
        end
        A(j,j)=A(j,j)-A(j,s:j-1)*v(s:j-1);
        A(j+1:u,j)=(A(j+1:u,j)-A(j+1:u,s:j-1)*v(s:j-1))/A(j,j);
    end
    for j = 1:n-1
        u = min(j+m,n);
        b(j+1:u)=b(j+1:u)-b(j)*A(j+1:u,j);
    end
    for i = 1:n
        b(i)=b(i)/A(i,i);
    end
    for j = n:-1:2
        s = max(j-m,1);
        b(s:j-1)=b(s:j-1)-b(j)*A(j,s:j-1)';
    end
end

function b = gauss(A,b)
n = length(A);
y = zeros(n,1);x = zeros(n,1);
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
for j = 1:n-1
    b(j+1:n)=b(j+1:n)-b(j)*A(j+1:n,j);
end
for j = n:-1:2
    b(j)=b(j)/A(j,j);
    b(1:j-1)=b(1:j-1)-b(j)*A(1:j-1,j);
end
b(1)=b(1)/A(1,1);
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

