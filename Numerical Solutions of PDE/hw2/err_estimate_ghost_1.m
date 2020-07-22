function err = err_estimate_ghost_1(U,real_sol,n,m,dt)
Err = zeros(n,m+1);
for k = 1:m+1
    for j = 1:n
    Err(j,k)=real_sol((j-0.5)/n,(k-1)*dt);
    end
end
Err=Err-U;
x =zeros(m*n+n,1);y =zeros(m*n+n,1);e =zeros(m*n+n,1);
for i = 1:n
    for j = 1:m+1
        x((j-1)*n+i)=(i-0.5)/n;
        y((j-1)*n+i)=(j-1)*dt;
        e((j-1)*n+i)=Err(i,j);
    end
end
scatter3(x,y,e);
Err=abs(Err);
err = max(max(Err));