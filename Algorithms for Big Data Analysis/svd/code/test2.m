function err_A = test2
m = 2048;
n = 512;
p = 20;
rand('seed',0);
A = randn(m,p)*randn(p,n);
[U, S, ~] = svd(A);
val = diag(S);
tic
p = zeros(1,n);
norm_A2 = norm(A,'fro')^2;
for i = 1:n
    p(i) = (norm(A(:,i)))^2/norm_A2;
end
[err1, s1, v1] = lineartime(A, 5, 200, p);
[err2, s2, v2] = lineartime(A, 10, 200, p);
[err3, s3, v3] = lineartime(A, 15, 200, p);
[err4, s4, v4] = lineartime(A, 20, 200, p);
toc


err_A = [err1; err2; err3; err4]

% err_Sigma = [norm(s1-val(1:5),'inf'); norm(s2-val(1:10),'inf'); norm(s3-val(1:15),'inf'); norm(s4-val(1:20),'inf')]
% err_vec = [norm(v1-U(:,1:5));norm(v2-U(:,1:10));norm(v3-U(:,1:15));norm(v4-U(:,1:20))]
fig = figure
hold off
subplot(2,2,1)
plot(1:5, val(1:5), 'o', 1:5, s1, '+', 'LineWidth', 1.2)
hold on
grid on
legend('builtin', 'lineartime')
title('r=5')
subplot(2,2,2)
plot(1:10, val(1:10), 'o', 1:10, s2, '+', 'LineWidth', 1.2)
legend('builtin', 'lineartime')
grid on
title('r=10')
subplot(2,2,3)
plot(1:15, val(1:15), 'o', 1:15, s3, '+', 'LineWidth', 1.2)
legend('builtin', 'lineartime')
grid on
title('r=15')
subplot(2,2,4)
plot(1:20, val(1:20), 'o', 1:20, s4,'+', 'LineWidth', 1.2)
legend('builtin', 'lineartime')
grid on
title('r=20')
hold off

print(fig, 'Report\test2.eps', '-depsc2'); 


fig = figure
hold off
[w1,ind] = sort(U(:,1));
w2 = v1(:,1);
w2 = w2(ind);
subplot(2,2,1)
plot(w1,'o','linewidth',2);
hold on
plot(w2,'x','linewidth',1);
legend('builtin','prototype')
title('r=5')

[w1,ind] = sort(U(:,1));
w2 = v2(:,1);
w2 = w2(ind);
subplot(2,2,2)
plot(w1,'o','linewidth',2);
hold on
plot(w2,'x','linewidth',1);
legend('builtin','prototype')
title('r=10')

[w1,ind] = sort(U(:,1));
w2 = v3(:,1);
w2 = w2(ind);
subplot(2,2,3)
plot(w1,'o','linewidth',2);
hold on
plot(w2,'x','linewidth',1);
legend('builtin','prototype')
title('r=15')

[w1,ind] = sort(U(:,1));
w2 = v4(:,1);
w2 = w2(ind);
subplot(2,2,4)
plot(w1,'o','linewidth',2);
hold on
plot(w2,'x','linewidth',1);
legend('builtin','prototype')
title('r=20')

print(fig, 'Report\test21.eps', '-depsc2'); 


fig = figure
subplot(2,2,1)
cor = zeros(5,1);
for i = 1:5
    cor(i) = abs(dot(U(:,i),v1(:,i))/norm(U(:,i))/norm(v1(:,i)));
end
xlabel('r');
ylabel('Correlation Coefficient')
plot(1:5,cor,'x','LineWidth', 1.2);
hold on
plot(cor);
grid on
title('r=5')
subplot(2,2,2)
cor = zeros(10,1);
for i = 1:10
    cor(i) = abs(dot(U(:,i),v2(:,i))/norm(U(:,i))/norm(v2(:,i)));
end
xlabel('r');
ylabel('Correlation Coefficient')
plot(1:10,cor,'x','LineWidth', 1.2);
hold on
plot(cor);
grid on
title('r=10')

subplot(2,2,3)
cor = zeros(15,1);
for i = 1:15
    cor(i) = abs(dot(U(:,i),v3(:,i))/norm(U(:,i))/norm(v3(:,i)));
end
xlabel('r');
ylabel('Correlation Coefficient')
plot(1:15,cor,'x','LineWidth', 1.2);
hold on
plot(cor);
grid on
title('r=15')


subplot(2,2,4)
cor = zeros(20,1);
for i = 1:20
    cor(i) = abs(dot(U(:,i),v4(:,i))/norm(U(:,i))/norm(v4(:,i)));
end
xlabel('r');
ylabel('Correlation Coefficient')
plot(1:20,cor,'x','LineWidth', 1.2);
hold on
plot(cor);
grid on
title('r=20')
hold off

print(fig, 'Report\test22.eps', '-depsc2'); 


end
