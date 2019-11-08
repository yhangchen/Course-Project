function [u,k] = diffusion_mixed_ghost11(a,theta,p,h,alpha,beta,q,r,n,dt,t_max)
%a,p,h,alpha,beta,q,r分别同报告里的定义。dx=1/n，dt分别为空间和时间步长。
mu = dt*n^2;
dig=zeros(n,1);off_dig=zeros(n,1);u = zeros(n,1);add=zeros(n,1);
%add为非齐次和边界值增加的项。dig和off_dig用来构建迭代矩阵的主对角和副对角。
for j = 1:n
    u(j)=h((j-0.5)/n);
end
k = 1;
t1 = clock;t2 = t1;
while etime(t2,t1) < t_max
for j = 1:n
    dig(j)=2*theta*mu*a((j-0.5)/n,(k-1)*dt)+1;
    off_dig(j)=-mu*a((j-0.5)/n,(k-1)*dt)*theta;
end
A = sparse([1:n,1:n-1,2:n],[1:n,2:n,1:n-1],[dig',off_dig(1:n-1)',...
    off_dig(2:n)'],n,n,3*n);
alpha1k = (2-alpha(k*dt)/n)/(2+alpha(k*dt)/n);
beta1k = (2-beta(k*dt)/n)/(2+beta(k*dt)/n);
alpha2k = 2/(2+alpha(k*dt)/n);beta2k = 2/(2+beta(k*dt)/n);
alpha1k1 = (2-alpha((k-1)*dt)/n)/(2+alpha((k-1)*dt)/n);
beta1k1 = (2-beta((k-1)*dt)/n)/(2+beta((k-1)*dt)/n);
alpha2k1 = 2/(2+alpha((k-1)*dt)/n);beta2k1 = 2/(2+beta((k-1)*dt)/n);
A(1,1)=A(1,1)-mu*a(0.5/n,(k-1)*dt)*theta*alpha1k;
A(n,n)=A(n,n)-mu*a(1-0.5/n,(k-1)*dt)*theta*beta1k;
for j=1:n
    dig(j)=-2*(1-theta)*mu*a((j-0.5)/n,(k-1)*dt)+1;
    off_dig(j)=mu*a((j-0.5)/n,(k-1)*dt)*(1-theta);
end
B = sparse([1:n,1:n-1,2:n],[1:n,2:n,1:n-1],[dig',off_dig(1:n-1)',...
    off_dig(2:n)'],n,n,3*n); 
B(1,1)=B(1,1)+mu*a(0.5/n,(k-1)*dt)*(1-theta)*alpha1k1;
B(n,n)=B(n,n)+mu*a(1-0.5/n,(k-1)*dt)*(1-theta)*beta1k1;
for j=1:n
    add(j)=dt*p((j-0.5)/n,(k-1)*dt);
end
add(1)=add(1)+mu*a(0.5/n,(k-1)*dt)*((1-theta)*alpha2k1*q((k-1)*dt)/n + ...
theta*alpha2k*q(k*dt)/n);
add(n)=add(n)+mu*a(1-0.5/n,(k-1)*dt)*((1-theta)*beta2k1*r((k-1)*dt)/n +...
theta*beta2k*r(k*dt)/n);
b = B*u+add;
u = sor(A,b,1.1);
k = k+1;t2 = clock;
end
k = k-1;
end



