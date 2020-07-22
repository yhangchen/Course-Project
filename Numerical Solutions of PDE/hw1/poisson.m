function u = poisson(f,g,n)
%n为每一维度的网格数。f，g定义与报告中相同。
%构建迭代矩阵。
A = sparse((n-1)^2,(n-1)^2);
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

%加入边值条件。
b(1) = b(1)+g(1/n,0)+g(0,1/n);
b(n-1)=b(n-1)+g(1/n,1)+g(0,1-1/n);
b((n-1)*(n-2)+1)=b((n-1)*(n-2)+1)+g(1,1/n)+g(1-1/n,0);
b((n-1)^2)=b((n-1)^2)+g(1,1-1/n)+g(1-1/n,1);
for i = 2:n-2
    b(i)=b(i)+g(0,i/n);
    b((i-1)*(n-1)+1)=b((i-1)*(n-1)+1)+g(i/n,0);
    b(i*(n-1))=b(i*(n-1))+g(i/n,1);
    b((n-1)*(n-2)+i)=b((n-1)*(n-2)+i)+g(1,i/n);
end

% 使用共轭梯度法求解Au=b;
x0 = b;
r0 = -A*b+b;
init_err = norm(r0);
err = init_err;
p0 = r0;
a0 = r0'*r0/(p0'*A*p0);
x0 = x0 + a0*p0;
r1 = r0 - a0*A*p0;
while err > 10^(-6)*init_err
    b0 = (r1'*r1)/(r0'*r0);
    p0 = r1 + b0*p0;
    p1 = A*p0;
    a0 = (r1'*r1)/(p0'*p1);
    x0 = x0 + a0*p0;
    r0 = r1;
    r1 = r0 - a0*p1;
    err = norm(r1);
end
u = x0;
%可视化解。调用此函数时不使用。
% x = zeros((n-1)^2,1);
% y = zeros((n-1)^2,1);
% for i = 1:n-1
% for j =1:n-1
% x((i-1)*(n-1)+j)=i/n;
% y((i-1)*(n-1)+j)=j/n;
% end
% end
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
% for i = 1:4*n
%     u_bd(i)=g(x_bd(i),y_bd(i));
% end
% x = [x' x_bd']';y = [y' y_bd']';u = [u' u_bd']';
% scatter3(x,y,u);
end
