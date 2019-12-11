function err = error_estimate(m,n,verbose,gridsize)
%input m: size of temporal grid, n: size of spatial grid, verbose: name of
%the method. gridsize: for numerical integration.
%output err: L2 error to the real solution.
    function z = real(x,y)
        z = exp(-2*pi^2)*sin(pi*x)*sin(pi*y);
    end

    function z = estimate(x,y,U0,n)
        x0 = floor(x*n);y0 = floor(y*n);
        if x0 == n
            x0 = n-1;
        end
        if y0 == n
            y0 = n-1;
        end
        xi = 2*n*(x-(x0+1/2)/n);mu = 2*n*(y-(y0+1/2)/n);
        z = U0(x0+1,y0+1)*(1-xi)*(1-mu)/4+U0(x0+1,y0+2)*(1-xi)*(1+mu)/4+...
        U0(x0+2,y0+1)*(1+xi)*(1-mu)/4+U0(x0+2,y0+2)*(1+xi)*(1+mu)/4;
    end
        
if strcmp(verbose,'explicit')
    [u,~] = explicit(m,n);
end
if strcmp(verbose,'implicit')
    u = timing(m,n,'gradient');
end
if strcmp(verbose,'crank')
    u = crank_improve(m,n);
end

U = zeros(n+1,n+1);
for i = 1:n-1
    for j = 1:n-1
        U(i+1,j+1) = u((i-1)*(n-1)+j);
    end
end
err = 0;
for i = 1:gridsize-1
    for j = 1:gridsize-1
        x = i/gridsize; y = j/gridsize;
        err = err + (real(x,y)-estimate(x,y,U,n))^2;
    end
end
err = sqrt(err)/gridsize;

% x = 0:1/100:1;
% y = 0:1/100:1;
% [xx,yy]=meshgrid(x,y);
% for i = 1:101
%     for j = 1:101
%         z(i,j) = real(x(i),y(j))-estimate(x(i),y(j),U,n);
%     end
% end
% plot3(xx,yy,z)
end