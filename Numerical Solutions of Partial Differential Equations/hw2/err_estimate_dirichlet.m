function err = err_estimate_dirichlet(U,real_sol,n,m,dt)
Err = zeros(n-1,m);
for k = 1:m
    for j = 1:n-1
    Err(j,k)=real_sol(j/n,k*dt);
    end
end
Err=Err-U(:,2:m+1);
x =zeros(m*(n-1),1);y =zeros(m*(n-1),1);e =zeros(m*(n-1),1);
for i = 1:n-1
    for j = 1:m
        x((j-1)*(n-1)+i)=i/n;
        y((j-1)*(n-1)+i)=j*dt;
        e((j-1)*(n-1)+i)=Err(i,j);
    end
end
scatter3(x,y,e);

%考虑l2收敛，则换为
% Err=Err(:,m);
% err = norm(Err);

%考虑l_infty收敛
Err=abs(Err);
err = max(max(Err));

end