n2 = 100;rand('seed',0);A2 = zeros(n2);b2=rand(n2,1);
for i = 1:n2-1
    A2(i,i)=10;A2(i,i+1)=1;A2(i+1,i)=1;
end
A2(n2,n2)=10;real2=A2\b2;
err21 = norm(improved_cho(A2,b2)-real2)
r21 = norm(A2*improved_cho(A2,b2)-b2)
err22 = norm(cho(A2,b2)-real2)
r22 = norm(cho(A2,b2)-b2)



n3 = 40;A3 = zeros(n3);b3=zeros(n3,1);real3=zeros(n3,1);
for i = 1:n3
    for j = 1:n3
        A3(i,j)=1/(i+j-1);
        b3(i)=b3(i)+1/(i+j-1);
    end
    real3(i)=1;
end
err31 = norm(improved_cho(A3,b3)-real3)
r31 = norm(A3*improved_cho(A3,b3)-b3)
err32 = norm(cho(A3,b3)-real3)
r32 = norm(A3*cho(A3,b3)-b3)


