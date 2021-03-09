function [u,k] = diffusion_mixed_ghost21(a,theta,p,h,alpha,beta,q,r,n,dt,t_max)
%a,p,h,alpha,beta,q,r分别同报告里的定义。dx=1/n，dt分别为空间和时间步长。
mu = dt*n^2;
dig=zeros(n+1,1);off_dig=zeros(n+1,1);u = zeros(n+1,1);add=zeros(n+1,1);
%add为非齐次和边界值增加的项。dig和off_dig用来构建迭代矩阵的主对角和副对角。
for j = 0:n
    u(j+1)=h(j/n);
end
t1 = clock; t2 = t1; k = 1;
while etime(t2,t1)<t_max
for j = 0:n
    dig(j+1)=2*theta*mu*a(j/n,(k-1)*dt)+1;
    off_dig(j+1)=-mu*a(j/n,(k-1)*dt)*theta;
end
A = sparse([1:n+1,1:n,2:n+1],[1:n+1,2:n+1,1:n],[dig',off_dig(1:n)',...
    off_dig(2:n+1)'],n+1,n+1,3*n+3);
A(1,1)=A(1,1)+mu*a(0,(k-1)*dt)*theta*2*alpha(k*dt)/n;
A(1,2)=A(1,2)-mu*a(0,(k-1)*dt)*theta;
A(n+1,n+1)=A(n+1,n+1)+mu*a(1,(k-1)*dt)*theta*2*beta(k*dt)/n;
A(n+1,n)=A(n+1,n)-mu*a(1,(k-1)*dt)*theta;
for j=0:n
    dig(j+1)=-2*(1-theta)*mu*a(j/n,(k-1)*dt)+1;
    off_dig(j+1)=mu*a(j/n,(k-1)*dt)*(1-theta);
end
B = sparse([1:n+1,1:n,2:n+1],[1:n+1,2:n+1,1:n],[dig',off_dig(1:n)',off_dig(2:n+1)'],...
n+1,n+1,3*n+3);
B(1,1)=B(1,1)-mu*a(0,(k-1)*dt)*(1-theta)*2*alpha((k-1)*dt)/n;
B(1,2)=B(1,2)+mu*a(0,(k-1)*dt)*(1-theta);
B(n+1,n+1)=B(n+1,n+1)-mu*a(1,(k-1)*dt)*(1-theta)*2*beta((k-1)*dt)/n;
B(n+1,n)=B(n+1,n)+mu*a(1,(k-1)*dt)*(1-theta);
for j=0:n
    add(j+1)=dt*p(j/n,(k-1)*dt);
end
add(1)=add(1)+mu*a(0,(k-1)*dt)*((1-theta)*2*q((k-1)*dt)/n+theta*2*q(k*dt)/n);
add(n+1)=add(n+1)+mu*a(1,(k-1)*dt)*((1-theta)*2*r((k-1)*dt)/n+theta*2*r(k*dt)/n);
b = B*u+add;
% u=A\b;
u = sor(A,b,1.1);
k = k + 1; t2 = clock;
end
k = k - 1;

end
