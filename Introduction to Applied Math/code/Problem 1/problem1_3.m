function X = problem1_3(threshold)
%input threshold for self-adaption
%output new grid
    function y = f(x)
        y = exp(-100*(x-0.5)^2);
    end
n = 4;
X = 0:1/4:1;
string = {};%for legend in the graph.
while n < threshold
    X = self_adaption(@(x) exp(-100*(x-0.5)^2), X, 0.5);
    n = length(X)-1;
    for i = 1:n
        h(i) = X(i+1)-X(i);
    end
    dig(1) = 1/h(1)+10^6;
    for i = 2:n
        dig(i) = 1/h(i-1)+1/h(i);
    end
    dig(n+1) = 1/h(n);
    for i = 1:n
        subdig(i) = -1/h(i);
    end
    x_axis = [1:n+1,1:n,2:n+1];
    y_axis = [1:n+1,2:n+1,1:n];
    A = sparse(x_axis,y_axis,[dig,subdig,subdig]);
    b(1) = f(X(1))*h(1)/2;
    for i = 1:n-1
        b(i+1) = f(X(i+1))*(h(i)+h(i+1))/2;
    end
    b(n+1) = f(X(n))*h(n)/2;
    for i = 2:n+1
        A(i,i) = A(i,i)-A(i-1,i)*A(i,i-1)/A(i-1,i-1);
        b(i) = b(i)-b(i-1)*A(i,i-1)/A(i-1,i-1);
        A(i,i-1) = 0;
    end
    mu(n+1) = b(n+1)/A(n+1,n+1);
    for i = n:-1:1
        mu(i) = (b(i)-A(i,i+1)*mu(i+1))/A(i,i);
    end
    scatter(X,mu)
    hold on
    string = [string,num2str(n)];
end
legend(string)
hold off
end


    