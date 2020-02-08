for i = 5:8
for j = 0:2:8
[A1(i-4,j/2+1),B1(i-4,j/2+1)]=poisson_choice(@f,10^(-j),2^i,'G-S');
[A2(i-4,j/2+1),B2(i-4,j/2+1)]=poisson_choice(@f,10^(-j),2^i,'Jacobi');
[A3(i-4,j/2+1),B3(i-4,j/2+1)]=poisson_choice(@f,10^(-j),2^i,'SOR');
end
end