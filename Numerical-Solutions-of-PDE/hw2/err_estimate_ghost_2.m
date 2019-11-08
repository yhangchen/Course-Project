function err = err_estimate_ghost_2(U,real_sol,n,m,dt)
Err = zeros(n+1,m);
for k = 1:m
    for j = 1:n+1
    Err(j,k)=real_sol((j-1)/n,k*dt);
    end
end
Err=Err-U(:,2:m+1);
x =zeros(m*(n+1),1);y =zeros(m*(n+1),1);e =zeros(m*(n+1),1);
for i = 1:n+1
    for j = 1:m
        x((j-1)*(n+1)+i)=(i-1)/n;
        y((j-1)*(n+1)+i)=j*dt;
        e((j-1)*(n+1)+i)=Err(i,j);
    end
end
scatter3(x,y,e);
Err=abs(Err);
err = max(max(Err));
end