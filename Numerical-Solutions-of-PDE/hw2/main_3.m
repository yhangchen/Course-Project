theta=0.8;n=40;dt=1/5000; t_max = 1;%调整此处
%不对m限制 

% 直接差分法用如下代码。这里left，right是最左边和最右边边界上值。
% [u,k,left,right] = diffusion_mixed_direct(@a,theta,@p,@h,@alpha,@beta,@q,@r,n,inf,dt,t_max);
% true = zeros(n-1,1);
% for i = 1:n-1
%     true(i)=real_sol(i/n,k*dt);
% end
% k
% err = max(norm(true-u,inf),max(abs(left-real_sol(0,k*dt)),abs(right-real_sol(1,k*dt))))

%第一种虚拟格点用如下代码
% [u,k]=diffusion_mixed_ghost11(@a,theta,@p,@h,@alpha,@beta,@q,@r,n,dt,t_max);
% true = zeros(n,1);
% for i = 1:n
%     true(i)=real_sol((i-1/2)/n,k*dt);
% end
% k
% err = norm(true-u,inf)

%第二种虚拟节点用以下代码
[u,k] = diffusion_mixed_ghost21(@a,theta,@p,@h,@alpha,@beta,@q,@r,n,dt,t_max);
true = zeros(n+1,1);
for i = 0:n
    true(i+1) = real_sol(i/n,k*dt);
end
k
err = norm(true-u,inf)
