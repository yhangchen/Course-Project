n2 = 1000;rand('seed',0);A2 = zeros(n2);b2=rand(n2,1);
for i = 1:n2-1
    A2(i,i)=10;A2(i,i+1)=1;A2(i+1,i)=1;
end
A2(n2,n2)=10;real2=A2\b2;
tic
err22 = norm(gauss_band(A2,b2,1)-real2,inf)
toc

tic
err2 = norm(least_square(A2,b2)-real2,inf)
toc

tic
err21 = norm(improved_cho(A2,b2)-real2,inf)
toc

