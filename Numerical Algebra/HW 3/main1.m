n1 = 84;A1 = zeros(n1);b1=zeros(n1,1);real1=zeros(n1,1);
for i = 1:n1-1
    A1(i,i)=6;A1(i,i+1)=1;A1(i+1,i)=8;
    if i>1
        b1(i)=15;
    end
    real1(i)=1;
end
A1(n1,n1)=6;b1(1)=7;b1(n1)=14;real1(n1)=1;
err1 = norm(least_square(A1,b1)-real1)
r1 = norm(A1*least_square(A1,b1)-b1)

rand('seed',0);
n2 = 100;A2 = zeros(n2);b2=rand(n2,1);
for i = 1:n2-1
    A2(i,i)=10;A2(i,i+1)=1;A2(i+1,i)=1;
end
A2(n2,n2)=10;real2=A2\b2;
err2 = norm(least_square(A2,b2)-real2)
r2 = norm(A2*least_square(A2,b2)-b2)


n3 = 40;A3 = zeros(n3);b3=zeros(n3,1);real3=zeros(n3,1);
for i = 1:n3
    for j = 1:n3
        A3(i,j)=1/(i+j-1);
        b3(i)=b3(i)+1/(i+j-1);
    end
    real3(i)=1;
end
err3 = norm(least_square(A3,b3)-real3)
r3 = norm(A3*least_square(A3,b3)-b3)


