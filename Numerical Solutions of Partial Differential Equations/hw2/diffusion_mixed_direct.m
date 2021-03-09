function [u,k,left,right] = diffusion_mixed_direct(a,theta,p,h,alpha,beta,q,r,n,m,dt,t_max)
%a,p,h,alpha,beta,q,r分别同报告里的定义。dx=1/n，dt分别为空间和时间步长。
%如果没有时间限制，将输出改为U，并取消下面与U有关注释。
mu = dt*n^2;
dig=zeros(n-1,1);off_dig=zeros(n-1,1);u = zeros(n-1,1);add=zeros(n-1,1);
%add为非齐次和边界值增加的项。dig和off_dig用来构建迭代矩阵的主对角和副对角。
for j = 1:n-1
    u(j)=h(j/n);
end

% U = zeros(n-1,m+1);
%U(:,1)=u;
% 做第三题时注释此行

t1 = clock;
t2 = t1;
k = 1;
left = h(0); right = h(1);
%记录最左侧和最右侧值。
while k <= m & etime(t2,t1) < t_max
for j = 1:n-1
    dig(j)=2*theta*mu*a(j/n,(k-1)*dt)+1;
    off_dig(j)=-mu*a(j/n,(k-1)*dt)*theta;
end
A = sparse([1:n-1,1:n-2,2:n-1],[1:n-1,2:n-1,1:n-2],[dig',off_dig(1:n-2)',...
    off_dig(2:n-1)'],n-1,n-1,3*n);
A(1,1)=A(1,1)-mu*a(1/n,(k-1)*dt)*theta/(1+alpha((k)*dt)/n);
A(n-1,n-1)=A(n-1,n-1)-mu*a(1-1/n,(k-1)*dt)*theta/(1+beta((k)*dt)/n);
for j=1:n-1
    dig(j)=-2*(1-theta)*mu*a(j/n,(k-1)*dt)+1;
    off_dig(j)=mu*a(j/n,(k-1)*dt)*(1-theta);
end
B = sparse([1:n-1,1:n-2,2:n-1],[1:n-1,2:n-1,1:n-2],[dig',off_dig(1:n-2)',...
    off_dig(2:n-1)'],n-1,n-1,3*n); 
B(1,1)=B(1,1)+mu*a(1/n,(k-1)*dt)*(1-theta)/(1+alpha((k-1)*dt)/n);
B(n-1,n-1)=B(n-1,n-1)+mu*a(1-1/n,(k-1)*dt)*(1-theta)/(1+beta((k-1)*dt)/n);
for j=1:n-1
    add(j)=dt*p(j/n,(k-1)*dt);
end
add(1)=add(1)+mu*a(1/n,(k-1)*dt)*(1-theta)/(1+alpha((k-1)*dt)/n)*q((k-1)*dt)/n + ...
mu*a(1/n,(k-1)*dt)*theta/(1+alpha((k)*dt)/n)*q(k*dt)/n;
add(n-1)=add(n-1)+mu*a(1-1/n,(k-1)*dt)*(1-theta)/(1+beta((k-1)*dt)/n)*r((k-1)*dt)/n +...
    mu*a(1-1/n,(k-1)*dt)*theta/(1+beta(k*dt)/n)*r(k*dt)/n;
b = B*u+add;
u = sor(A,b,1.1);
left=(u(1)+q(k*dt)/n)/(1+alpha(k*dt)/n);
right=(u(n-1)+r(k*dt)/n)/(1+beta(k*dt)/n);
% U(:,k+1)=u;
k = k+1; t2 = clock;
end
k = k-1;
end



