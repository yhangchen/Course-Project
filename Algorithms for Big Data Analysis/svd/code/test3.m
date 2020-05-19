function test3(t,m,n,k)
[A, U, val] = PCAtestmatrix(m, n, t);
fprintf('prototype')
tic
[err1, s1, u1] = prototype(A, k, 0);
err1
toc
fprintf('lineartime')
tic
p = zeros(1,n);
norm_A2 = norm(A,'fro')^2;
for i = 1:n
    p(i) = (norm(A(:,i)))^2/norm_A2;
end
[err2, s2, u2] = lineartime(A, k, 10*k, p);
err2
toc

fig = figure;
hold off
subplot(2,2,1)
semilogy(1:k, val(1:k), 'o', 1:k, s1, 'x', 1:k, s2, 'x','LineWidth', 1.2)
hold on
grid on
xlabel('$k$-th singular value','interpreter','latex')
legend('groundtruth', 'prototype','lineartime','interpreter','latex')
ylabel('singular value','interpreter','latex')
subplot(2,2,2)
plot(1:k, abs(-val(1:k)+s1), 'x', 1:k, abs(-val(1:k)+s2), 'x','LineWidth', 1.2)
grid on
xlabel('$k$-th singular value','interpreter','latex')
legend('prototype','lineartime','interpreter','latex')
ylabel('absolute error','interpreter','latex')
hold off





cor1 = zeros(k,1);
cor2 = zeros(k,1);
for i = 1:k
    cor1(i) = abs(dot(U(:,i),u1(:,i))/norm(U(:,i))/norm(u1(:,i)));
    cor2(i) = abs(dot(U(:,i),u2(:,i))/norm(U(:,i))/norm(u2(:,i)));
end

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

e21 = U(:,1)-u2(:,1);
e22 = U(:,1)+u2(:,1);
if norm(e21) < norm(e22)
    e2 = e21;
else
    e2 = e22;
end

plot(1:m,e1,'x',1:m,e2,'x');
grid on
xlabel('index of the vector','interpreter','latex')
ylabel('error of the first singular vector','interpreter','latex')
legend('prototype','lineartime','interpreter','latex')

subplot(2,2,4)
plot(1:k,cor1,'-x',1:k,cor2,'-x','LineWidth', 1.2);
xlabel('$k$-th singular vector','interpreter','latex')
legend('prototype','lineartime','interpreter','latex')
ylabel('correlation','interpreter','latex')
hold on
grid on
hold off

set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','points');
set(gcf,'PaperPosition', [0 0 640 400]);
set(gcf,'Position', [0 0 640 400]);



name = 'Report\'+string(t)+'_'+string(k)+'.eps';
print(fig, name, '-depsc2'); 

end

