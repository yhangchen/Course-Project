clear all
% comparison of two criteria, model 1/2
n = 30; model = 1;
if model == 1
    S0 = model1(n);
else
    [S0,flag] = model2(n);
    if flag ~= 1
        error('model 2 has no solution, please retry');
    end
end
test = mvnrnd(zeros(n,1),S0,100);
training = mvnrnd(zeros(n,1),S0,100);
S_train = cov(training); 
S_test = cov(test);
rhos = 0.001:0.005:0.2;
for j = 1:length(rhos)
    rho = rhos(j)
    [X1, ~] = spinv_glasso(S_train, rho);
    [X2, ~] = spinv_admm(S_train, 1/rho);
%     l1_norm(j) = sum(sum(abs(X1)))-sum(sum(abs(X2)));
%     F_norm(j) = norm(X1,'fro')-norm(X2,'fro');
%     op_norm(j) = norm(X1)-norm(X2);
%     zero_norm(j) = sum(sum(abs(X1)>eps))-sum(sum(abs(X2)>eps));
%     val1(j) = trace(S_test'*X1)-log(det(X1));
%     val2(j) = trace(S_test'*X2)-log(det(X2));
val1(j) = norm(S_test'*X1-eye(n),'fro');
val2(j) = norm(S_test'*X2-eye(n),'fro');
end
hold off
hold on
plot(rhos,val1,'-x','linewidth',2)
plot(rhos,val2,'-x','linewidth',2)
xlabel('$\rho=1/\sigma$','interpreter','latex','Fontsize',12);
ylabel('F-norm loss','Fontsize',12);
legend('Criterion 1', 'Criterion 2','Fontsize',12);

