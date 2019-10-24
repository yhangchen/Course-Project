function err = err_estimate_r(f,g,h,m,n)
%h为真解。
u = poisson_radius(f,g,m,n);
v = zeros((n-1)*m,1);
for i = 1:n-1
    for j = 1:m
        v((i-1)*m+j)=h((j-0.5)/m,pi*i/2/n);
    end
end
r = u-v;
err = norm(r,inf);
x = zeros((n-1)*m,1);y=zeros((n-1)*m,1);
for j = 1:m
    for k = 1:n-1
        x((k-1)*m+j)=(j-0.5)/m*cos(pi*k/2/n);
        y((k-1)*m+j)=(j-0.5)/m*sin(pi*k/2/n);
    end
end

% scatter3(x,y,r);
%单独使用此程序绘图时取消注释。
end