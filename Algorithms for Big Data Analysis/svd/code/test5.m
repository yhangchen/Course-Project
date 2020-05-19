function test5(k)
argin = load('matlab.mat');
A = argin.A;
m = 1e4; n=1e4;

S0 = zeros(n,1);
for i=1:20
    S0(i)=10^(-4*(i-1)/19);
end

for i=21:n
   S0(i)=(10^(-4))/(i-20)^(1/10);
end


fprintf('svds')
tic
[U,S,~] = svds(A,k);
toc
val = diag(S);
fprintf('prototype')
tic
[err1, s1, u1] = prototype(A, k, 0);
toc
err1

fig = figure;
hold off
subplot(2,2,1)
semilogy(1:k, S0(1:k), 'o', 1:k, val, 'x', 1:k, s1, 'x','LineWidth', 1.2)
hold on
grid on
xlabel('$k$-th singular value','interpreter','latex')
legend('groundtruth','svds', 'prototype' ,'interpreter','latex')
ylabel('singular value','interpreter','latex')
subplot(2,2,2)
plot(1:k, abs(-S0(1:k)+val), 'x',1:k, abs(-S0(1:k)+s1), 'x','LineWidth', 1.2)
grid on
xlabel('$k$-th singular value','interpreter','latex')
legend('svds','prototype' ,'interpreter','latex')
ylabel('absolute error','interpreter','latex')
hold off





cor1 = zeros(k,1);
% cor2 = zeros(k,1);
for i = 1:k
    cor1(i) = abs(dot(U(:,i),u1(:,i))/norm(U(:,i))/norm(u1(:,i)));
%     cor2(i) = abs(dot(U(:,i),u2(:,i))/norm(U(:,i))/norm(u2(:,i)));
end


% in the following, we regard svds as groundtruth

subplot(2,2,3)

% since u and -u are both singular vectors.
% we need to regard u and -u identical.
e11 = U(:,1)-u1(:,1);
e12 = U(:,1)+u1(:,1);
if norm(e11) < norm(e12)
    e1 = e11;
else
    e1 = e12;
end

plot(1:m,e1,'x');
grid on
xlabel('index of the vector','interpreter','latex')
ylabel('error of the first singular vector','interpreter','latex')
legend('prototype' ,'interpreter','latex')

subplot(2,2,4)
plot(1:k,cor1,'-x','LineWidth', 1.2);
xlabel('$k$-th singular vector','interpreter','latex')
legend('prototype' ,'interpreter','latex')
ylabel('correlation','interpreter','latex')
hold on
grid on
hold off

set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','points');
set(gcf,'PaperPosition', [0 0 640 400]);
set(gcf,'Position', [0 0 640 400]);



name = 'Report\large'+string(k)+'.eps';
print(fig, name, '-depsc2'); 

end

