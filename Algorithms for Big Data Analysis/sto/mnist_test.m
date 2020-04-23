function mnist_test(batchsize, lambda)
%% params setting
max_iter = 30;
min_iter = 10;
% batchsize = 200;
stopping = 0;
lambda0 = lambda/28/28; % resize lambda
% The parameters are the same as Pytorch's default setting for 
% torch.optim.Adagrad and torch.optim.Adam
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta1 = 0.9;
beta2 = 0.999;
eps_adam = 1e-8;
lr_adam = 0.001;

eps_adagrad = 1e-10;
lr_adagrad = 0.001;

m_sgd = 60000/batchsize;
beta_sgd = 10/m_sgd;
lr_sgd_0 = 1e-4; lr_sgd_1 = lr_sgd_0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
data = mnist_loader;
fprintf('Data loaded')

x = data.train_img;
x = reshape(x, [28*28, 60000])/255;
y = data.train_lbl;
y = mod(y,2) == 0;
y = 2*y - 1;
y = y';

%% AdaGrad
w00 = randn([784, 1])/batchsize;

w1 = w00; w2 = w00; w3 = w00;
r = 0;
t1 = 0;
[value1, grad1, loss1] = objective(x, y, w1, lambda0);
value0 = 0;
Value1 = [value1]; Loss1 = [loss1]; Grad1 = [norm(grad1,inf)];
FValue1 = [value1]; FLoss1 = [loss1]; FGrad1 = [norm(grad1,inf)];
while t1 < max_iter
    if t1 > min_iter
        if abs(value0-value1) < stopping; break; end
    end
    value0 = value1;
    t1 = t1 + 1;
    value11 = 0; loss11 = 0; grad11 = 0;
    for i = 1:60000/batchsize
    ind = randi(60000, [batchsize ,1]);
    x0 = x(:,ind);
    y0 = y(:,ind);
    [value1, grad1, loss1, w1, r] = AdaGrad(@objective, x0, y0, w1, lambda0, r, t1, lr_adagrad, eps_adagrad);
    value11 = value11 + value1;
    grad11 = grad11 + norm(grad1,inf);
    loss11 = loss11 + loss1;
    end
    
%     [fvalue1, fgrad1, floss1] = objective(x, y, w1, lambda/60000);
    
    Value1 = [Value1 value11*batchsize/60000];
    Grad1 = [Grad1 grad11*batchsize/60000];
    Loss1 = [Loss1 loss11*batchsize/60000];
 
%     FValue1 = [FValue1 fvalue1];
%     FGrad1 = [FGrad1 norm(fgrad1,inf)];
%     FLoss1 = [FLoss1 floss1];
% 
end


%% ADAM
% w2 = randn([784, 1])/batchsize;
m = 0;
v = 0;
t2 = 0;

[value1, grad1, loss1] = objective(x, y, w2, lambda0);
value0 = 0;
Value2 = [value1]; Loss2 = [loss1]; Grad2 = [norm(grad1,inf)];
while t2 < max_iter
    if t2 > min_iter
        if abs(value0-value1) < stopping; break; end
    end
    value0 = value1;
    t2 = t2 + 1;
    value11 = 0; loss11 = 0; grad11=0;
    for i = 1:60000/batchsize
    ind = randi(60000, [1,batchsize]);
    x0 = x(:,ind);
    y0 = y(:,ind);
    it = (t2 - 1)*60000/batchsize + i;
    [value1, grad1, loss1, w2, m, v] = ADAM(@objective, x0, y0, w2, lambda0, m, v, it, lr_adam, eps_adam, beta1, beta2);
    value11 = value11 + value1;
    grad11 = grad11 + norm(grad1,inf);
    loss11 = loss11 + loss1;
    end
    Value2 = [Value2 value11*batchsize/60000];
    Grad2 = [Grad2 grad11*batchsize/60000];
    Loss2 = [Loss2 loss11*batchsize/60000];
end

%% SGD-BB
[Value3, Grad3, Loss3, w3, t3] = SGD_BB(@objective, x, y, w3, lambda0, m_sgd, beta_sgd, lr_sgd_0, lr_sgd_1, max_iter, batchsize, 60000, stopping);

%% plot the figures
fig = figure(1);
s1 = 'AdaGrad iter='+string(t1); s2 = 'Adam iter='+ string(t2); s3 = 'SGD-BB iter='+ string(t3);
hold off
plot(Value1,'-x','linewidth',2)
hold on
grid on
plot(Value2,'-x','linewidth',2)
plot(Value3,'-x','linewidth',2)
xlabel('epoch')
ylabel('value of $f$','interpreter','latex')
legend(s1,s2,s3)
title('$\lambda='+string(lambda)+'$','interpreter','latex')
name = 'Report\mnist_'+string(lambda)+'_'+string(batchsize)+'.eps';
print(fig, name, '-depsc2'); 
hold off

fig = figure(2);
s1 = 'AdaGrad iter='+string(t1); s2 = 'Adam iter='+ string(t2); s3 = 'SGD-BB iter='+ string(t3);
hold off
plot(Loss1,'-x','linewidth',2)
hold on
grid on
plot(Loss2,'-x','linewidth',2)
plot(Loss3,'-x','linewidth',2)
xlabel('epoch')
ylabel('0-1 loss','interpreter','latex')
legend(s1,s2,s3)
title('$\lambda='+string(lambda)+'$','interpreter','latex')
name = 'Report\loss_mnist_'+string(lambda)+'_'+string(batchsize)+'.eps';
print(fig, name, '-depsc2'); 
hold off

fig = figure(3);
s1 = 'AdaGrad iter='+string(t1); s2 = 'Adam iter='+ string(t2); s3 = 'SGD-BB iter='+ string(t3);
hold off
plot(Grad1,'-x','linewidth',2)
hold on
grid on
plot(Grad2,'-x','linewidth',2)
plot(Grad3,'-x','linewidth',2)
xlabel('epoch')
ylabel('$\ell_\infty$ norm of gradient','interpreter','latex')
legend(s1,s2,s3)
title('$\lambda='+string(lambda)+'$','interpreter','latex')
name = 'Report\norm_mnist_'+string(lambda)+'_'+string(batchsize)+'.eps';
print(fig, name, '-depsc2'); 
hold off


%% test
x0 = data.test_img;
x0 = reshape(x0, [28*28, 10000]);
y0 = data.test_lbl;
y0 = mod(y0,2) == 0;
y0 = 2*y0 - 1;
y0 = y0';

acc_adagrad = mean((y0.*(w1'*x0) > 0))
acc_adam = mean((y0.*(w2'*x0) > 0))
acc_sgdbb = mean((y0.*(w3'*x0) > 0))