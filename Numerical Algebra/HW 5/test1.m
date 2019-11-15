function [err, k, t] = test1(n)
A = zeros(n);b=zeros(n,1);real=zeros(n,1);
for i = 1:n
    for j = 1:n
        A(i,j)=1/(i+j-1);
        b(i)=b(i)+1/(i+j-1)/3;
    end
    real(i)=1/3;
end
[x,k,t] = conjugate_gradient(A, b, 1e-14, inf);
err = norm(x-real);

% use the following code to generate figure.
% x = 1:100;
% for i = 1:100
% [y(i),z(i),w(i)] = test1(i);
% end
% plot(x,y)
% plot(x,z)
% plot(x,w)