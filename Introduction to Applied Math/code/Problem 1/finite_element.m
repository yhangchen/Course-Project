function mu = finite_element(f,X)
%input d^2u/dx^2 = -f, and grid X, eg X = 0:0.1:1
%output value of u at each grid point.
n = length(X)-1;
for i = 1:n
    h(i) = X(i+1)-X(i);
end
for i = 1:n-1
    dig(i) = 1/h(i)+1/h(i+1);
end
subdig = zeros(1,n-2);
for i = 2:n-1
    subdig(i-1) = -1/h(i);
end
x_axis = [1:n-1,1:n-2,2:n-1];
y_axis = [1:n-1,2:n-1,1:n-2];
A = sparse(x_axis,y_axis,[dig,subdig,subdig]);
%solve linear equation.
for i = 1:n-1
    b(i) = f(X(i+1))*(h(i)+h(i+1))/2;
end
for i = 2:n-1
    A(i,i) = A(i,i)-A(i-1,i)*A(i,i-1)/A(i-1,i-1);
    b(i) = b(i)-b(i-1)*A(i,i-1)/A(i-1,i-1);
    A(i,i-1) = 0;
end
mu(n) = b(n-1)/A(n-1,n-1);
for i = n-2:-1:1
    mu(i+1) = (b(i)-A(i,i+1)*mu(i+2))/A(i,i);
end
%add boundary part the solution.
mu(1) = 0; mu(n+1) = 0;
end

    