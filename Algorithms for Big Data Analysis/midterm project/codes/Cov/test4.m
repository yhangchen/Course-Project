clear all
hold off
% find optimal sigma=1/rho; under model 1/2
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
ms = [20,30,40,80,160];
rel_diff = @(x) x/x(1);
S_test = cov(test);
rhos = 0.01:0.05:1;%0.01:0.05:1;
for j = 1:length(ms)
    m = ms(j);
    training = mvnrnd(zeros(n,1),S0,m);
    S = cov(training); 
    val{j} = zeros(size(rhos));
    for i = 1:length(rhos)
        rho = rhos(i);
        [X1, ~] = spinv_glasso(S, rho);
        val{j}(i) = trace(S_test'*X1)-log(det(X1));
    end
end
fig = figure;
hold on
subplot(1,2,1)
hold on
plot(rhos,val{1},'-x','linewidth',2)
plot(rhos,val{2},'-x','linewidth',2)
plot(rhos,val{3},'-x','linewidth',2)
plot(rhos,val{4},'-x','linewidth',2)
plot(rhos,val{5},'-x','linewidth',2)
xlabel('$1/\sigma$','interpreter','latex','Fontsize',12);
ylabel('Likelihood loss','Fontsize',12);
legend('likelihood loss, m=20','likelihood loss, m=30','likelihood loss, m=40',...
    'likelihood loss, m=80','likelihood loss, m=160','Fontsize',12)

subplot(1,2,2)
hold on
plot(rhos,rel_diff(val{1}),'-x','linewidth',2)
plot(rhos,rel_diff(val{2}),'-x','linewidth',2)
plot(rhos,rel_diff(val{3}),'-x','linewidth',2)
plot(rhos,rel_diff(val{4}),'-x','linewidth',2)
plot(rhos,rel_diff(val{5}),'-x','linewidth',2)
legend('relative likelihood loss, m=20','relative likelihood loss, m=30','relative likelihood loss, m=40',...
    'relative likelihood loss, m=80','relative likelihood loss, m=160','Fontsize',12)
xlabel('$1/\sigma$','interpreter','latex','Fontsize',12);
ylabel('Likelihood loss','Fontsize',12);
hold off
