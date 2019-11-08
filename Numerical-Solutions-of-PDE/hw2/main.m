%需要适当调整这里的内容以输出报告中结果。
theta=0.5;
for t = 6:2:16
n=100;m=t*100;dt=0.01;
tic
% U=diffusion_dirichlet(@a,theta,@p,@h,@f,@g,n,m,dt);
U = diffusion_mixed_ghost2(@a,theta,@p,@h,@alpha,@beta,@q,@r,n,m,dt);
toc
err_estimate_ghost_2(U,@real_sol,n,m,dt)
%换成err_estimate_ghost_2
end

