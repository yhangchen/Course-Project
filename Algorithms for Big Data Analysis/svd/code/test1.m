function test1
m = 2048;
n = 512;
p = 20;
q = 1;
randn('seed',0)
A1 = randn(m,p);
randn('seed',0)
A2 = randn(p,n);
A = A1 * A2;
[U, S, ~] = svd(A);
val = diag(S);
tic
[err1, s1, v1] = prototype(A, 5, q);
[err2, s2, v2] = prototype(A, 10, q);
[err3, s3, v3] = prototype(A, 15, q);
[err4, s4, v4] = prototype(A, 20, q);
toc


err_A = [err1; err2; err3; err4]

figure
hold off
subplot(2,2,1)
plot(1:5, val(1:5), 'o', 1:5, s1, '+', 'LineWidth', 1.2)
hold on
grid on
legend('builtin', 'prototype')
title('r=5')
subplot(2,2,2)
plot(1:10, val(1:10), 'o', 1:10, s2, '+', 'LineWidth', 1.2)
legend('builtin', 'prototype')
grid on
title('r=10')
subplot(2,2,3)
plot(1:15, val(1:15), 'o', 1:15, s3, '+', 'LineWidth', 1.2)
legend('builtin', 'prototype')
grid on
title('r=15')
subplot(2,2,4)
plot(1:20, val(1:20), 'o', 1:20, s4,'+', 'LineWidth', 1.2)
legend('builtin', 'prototype')
grid on
title('r=20')
hold off


%for the sake of clarity, we sort the vector's items.
%%
figure
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

%%

figure
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
ylim([0.91,1.01]);
grid on
title('r=20')



hold off
end
