function z = Enneper(x,y)
% returns the surface's value on the point (x,y)
F = @(X) [X(1) + X(1)*X(2)^2 - X(1)^3/3 - x; -X(2) - X(1)^2*X(2) + X(2)^3/3 - y];
U = fsolve(F,[0,0],optimset('FunValCheck', 'off', 'Display', 'off'));
z = U(1)^2-U(2)^2;
end