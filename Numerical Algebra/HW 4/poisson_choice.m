function [err,iter] = poisson_choice(f,e,n,choice)
%构建迭代矩阵。
A = spalloc((n-1)^2,(n-1)^2,4*n^2);
for i = 1:(n-1)^2
    A(i,i)=2*(1+e);
end
for i = 1:(n-1)*(n-2)
    A(i,i+n-1)=-1;A(i+n-1,i)=-1;
end
for i = 1:(n-1)^2-1
     A(i,i+1)=-e;A(i+1,i)=-e;
end
for i = 1:n-2
    A((n-1)*i,(n-1)*i+1)=0;A((n-1)*i+1,(n-1)*i)=0;
end
b = zeros((n-1)^2,1);
for i = 1:n-1
    for j = 1:n-1
        b((i-1)*(n-1)+j)=f(i/n,j/n,e)/n^2;
    end
end
iter = 0;
%以下我们均用内置\来解三角方程，因为稀疏矩阵索引较慢，如果写开，即使采用稀疏格式，也很难在合理时间内得到结论。
if strcmp(choice,'G-S')
    t0 = clock;
    L = tril(A);
    U = -triu(A,1);
    u = b;
    init_err = norm(A*u-b);
    err = init_err;
    while err > 10^(-6)*init_err
        u = L\(U*u + b);
        err = norm(A*u-b);
        iter = iter + 1;
    end
    t = etime(clock,t0);
end
if strcmp(choice,'Jacobi')
    t0 = clock;
    u = b;
    init_err = norm(A*u-b);
    err = init_err;
    while err > 10^(-6)*init_err
        u = u - (A * u - b)/2/(1 + e);
        err = norm(A*u-b);
        iter = iter + 1;
    end
    t = etime(clock,t0);
end
if strcmp(choice,'SOR')
    t0 = clock;
    D = 2*(1+e)*speye((n-1)^2);
    U = -triu(A,1);
    w = 1.8; 
    u = b;
    init_err = norm(A*u-b);
    err = init_err;
    while err > 10^(-6)*init_err
        u = (D-w*U')\(((1-w)*D+w*U)*u+w*b);
        err = norm(A*u-b);
        iter = iter + 1;
    end
    t = etime(clock,t0);
end
u = A\b;

%可视化解。
x = zeros((n-1)^2,1);
y = zeros((n-1)^2,1);
e = zeros((n-1)^2,1);
for i = 1:n-1
for j =1:n-1
x((i-1)*(n-1)+j)=i/n;
y((i-1)*(n-1)+j)=j/n;
e((i-1)*(n-1)+j)=u((i-1)*(n-1)+j)-real_sol(i/n,j/n);
end
end
err = norm(e,inf);
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




