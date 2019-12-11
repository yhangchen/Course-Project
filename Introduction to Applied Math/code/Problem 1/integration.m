function res = integration(f,X)
%input function f, and grid X, eg X = 0:0.1:1
%output vector for problem 1.2(e2).(not for problem 1.1)
n = length(X)-1;
for i = 1:n
    res(i) = (X(i+1)-X(i))*sqrt((f(X(i+1))^2+f(X(i))^2)*(X(i+1)-X(i))/2);
end
end
