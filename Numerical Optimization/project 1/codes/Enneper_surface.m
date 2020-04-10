function z = Enneper_surface(nx,ny,data)
%{
 Enneper surface with grid size nx x ny, if data is not specified
plot
 input:
     nx, ny
 output
     z:
     matrix of (nx+1)x(ny+1), value on the grid points, and plot the surface.
%}
xs = 1/nx*[0:nx] - 0.5;
ys = 1/ny*[0:ny] - 0.5;
[y,x] = meshgrid(ys,xs);
z = zeros(nx+1,ny+1);
if nargin < 3
    for i = 1:nx+1
        for j = 1:ny+1
            z(i,j) = Enneper(xs(i),ys(j));
        end
    end
else
    z(2:nx,2:ny) = data.in;
    z(1,:) = data.left;
    z(nx+1,:) = data.right;
    z(:,1) = data.bottom;
    z(:,ny+1) = data.top;
end
s = mesh(x,y,z);
s.FaceColor = 'flat';
s.EdgeColor = 'c';
end