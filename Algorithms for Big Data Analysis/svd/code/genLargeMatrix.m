% Generate large matrix data using DCT.
m=1E4;                                 % size of matrix
n=1E4;

% from https://github.com/WenjianYu/rSVD-single-pass/
% I have trouble loading dat file, it seems that a 1e4 x 1e4 matrix will
% become a 1e4 x 5e3 matrix.

% hence I set smaller m,n and use a dense
% matrix to store it.

A = zeros(m,n);
% % % % % type 1% % % % %
S = zeros(n,1);
for i=1:20
    S(i)=10^(-4*(i-1)/19);
end

for i=21:n
   S(i)=(10^(-4))/(i-20)^(1/10);
end
% % % % % type 2% % % % %
% for i = 1:n
%    S(i) = i^(-2);
% end

% % % % % type 3% % % % %
%for i = 1:n
 %   S(i) = i^(-3);
%end

s = spdiags(S,0,n,m);

for i=1:m
    E=zeros(m,1);
    E(i)=1;
    e=idct(E);
    f=s*e;
    F=idct(f);
    A(i,:) = f';
end
