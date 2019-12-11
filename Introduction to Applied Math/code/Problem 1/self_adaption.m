function new_X = self_adaption(f, X, alpha)
%input d^2u/dx^2 = -f, grid X, threshold alpha.
%output new grid
n = length(X)-1;
res = integration(f,X);
mx = max(res);
new_point = [];
for i = 1:n
    if res(i) > alpha*mx
         new_point = [new_point, (X(i)+X(i+1))/2];
         %add new point.
    end
end
new_X = sort([X,new_point]);