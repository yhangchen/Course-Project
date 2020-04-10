function [f,g] = test_minima_surface(v, nx, ny, bottom, top, left, right)
%{
\min \{\sum\left(
        f_{i,j}^L(v)+f_{i,j}^U(v)
    \right):v\in \mathbb{R}\}
where
        f_{i,j}^L(v)=\frac{h_xh_y}{2}\left[
            1+\left(\frac{v_{i+1,j}-v_{i,j}}{h_x}\right)^2+\left(\frac{v_{i,j+1}-v_{i,j}}{h_x}\right)^2
        \right]^{1/2}\\
        f_{i,j}^U(v)=\frac{h_xh_y}{2}\left[
            1+\left(\frac{v_{i-1,j}-v_{i,j}}{h_x}\right)^2+\left(\frac{v_{i,j-1}-v_{i,j}}{h_x}\right)^2
        \right]^{1/2}
    
computes the function and gradient of the minimal surface area problem.
    input
       nx
           the number of grid points in the first
       ny
           the number of grid points in the second
       v
           vector size (nx-1)x(ny-1), encode the 
       
       bottom 
           vector dimension nx + 1.
           bottom must contain boundary data beginning
                with the lower left corner of the domain.

       top
           vector dimension nx + 1.
             On entry top must contain boundary data beginning with
                the upper left corner of the domain.

       left
           vector dimension nx + 1.
             left must contain boundary data beginning with
                the lower left corner of the domain.

       right 
           vector dimension nx + 1.
             right must contain boundary data beginning with
                the lower right corner of the domain.

   output
       f
           function value

       g
           vector of dimension (nx-1)*(ny-1), encode the gradient.
        modified from Minpack 2 Fortran code.
%}
hx = 1/nx; hy = 1/ny;
g = zeros(nx-1,ny-1);

V = zeros(nx+1,ny+1);
V(2:nx,2:ny) = reshape(v,nx-1,ny-1);
V(:,1) = bottom; V(:,ny+1) = top; V(1,:) = left; V(nx+1,:) = right;
dvdx = (V(2:nx+1,:)-V(1:nx,:))/hx;
dvdy = (V(:,2:ny+1)-V(:,1:ny))/hy;
fl = sqrt(1+dvdx(1:nx,1:ny).^2+dvdy(1:nx,1:ny).^2);
fu = sqrt(1+dvdx(1:nx,2:ny+1).^2+dvdy(2:nx+1,1:ny).^2);
f = hx*hy*(sum(sum(fl+fu)))/2;
for i=2:nx
    for j=2:ny
    g(i-1,j-1)=(1/hx)*(1./fl(i-1,j)+1./fu(i-1,j-1)).*dvdx(i-1,j)...
    -(1/hx)*(1./fl(i,j)+1./fu(i,j-1)).*dvdx(i,j)...
    +(1/hy)*(1./fl(i,j-1)+1./fu(i-1,j-1)).*dvdy(i,j-1)...
    -(1/hy)*(1./fl(i,j)+1./fu(i-1,j)).*dvdy(i,j);
    end
end
g = hx*hy*g(:)/2;
end
