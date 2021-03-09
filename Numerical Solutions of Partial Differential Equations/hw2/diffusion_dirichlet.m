function U = diffusion_dirichlet(a,theta,p,h,f,g,n,m,dt)
%a,p,h,f,g分别同报告里的定义。dx=1/n，dt分别为空间和时间步长。
mu = dt*n^2;N=0;%计算迭代步数
U = zeros(n-1,m+1);
dig=zeros(n-1,1);off_dig=zeros(n-1,1);u = zeros(n-1,1);add=zeros(n-1,1);
%add为非齐次和边界值增加的项。dig和off_dig用来构建迭代矩阵的主对角和副对角。
for j = 1:n-1
    u(j)=h(j/n);
end
U(:,1)=u;
for k = 1:m
for j = 1:n-1
    dig(j)=2*theta*mu*a(j/n,(k-1)*dt)+1;
    off_dig(j)=-mu*a(j/n,(k-1)*dt)*theta;
end
A = sparse([1:n-1,1:n-2,2:n-1],[1:n-1,2:n-1,1:n-2],[dig',off_dig(1:n-2)',off_dig(2:n-1)'],...
n-1,n-1,3*n);
for j=1:n-1
    dig(j)=-2*(1-theta)*mu*a(j/n,(k-1)*dt)+1;
    off_dig(j)=mu*a(j/n,(k-1)*dt)*(1-theta);
end
B = sparse([1:n-1,1:n-2,2:n-1],[1:n-1,2:n-1,1:n-2],[dig',off_dig(1:n-2)',off_dig(2:n-1)'],...
n-1,n-1,3*n); 
for j=1:n-1
    add(j)=dt*p(j/n,(k-1)*dt);
end
add(1)=add(1)+mu*a(1/n,(k-1)*dt)*theta*f((k)*dt)+mu*a(1/n,(k-1)*dt)*(1-theta)*f((k-1)*dt);
add(n-1)=add(n-1)+mu*a(1-1/n,(k-1)*dt)*theta*g((k)*dt)+mu*a(1-1/n,(k-1)*dt)*(1-theta)*g((k-1)*dt);
b = B*u+add;
% 使用共轭梯度法求解Au=b;
x0 = b;
r0 = -A*b+b;
init_err = norm(r0);
err = init_err;
p0 = r0;
a0 = r0'*r0/(p0'*A*p0);
x0 = x0 + a0*p0;
r1 = r0 - a0*A*p0;
while err > 10^(-15)*init_err
    b0 = (r1'*r1)/(r0'*r0);
    p0 = r1 + b0*p0;
    p1 = A*p0;
    a0 = (r1'*r1)/(p0'*p1);
    x0 = x0 + a0*p0;
    r0 = r1;
    r1 = r0 - a0*p1;
    err = norm(r1);
    N = N+1;
end
u = x0;
U(:,k+1)=u;
end
N
end



