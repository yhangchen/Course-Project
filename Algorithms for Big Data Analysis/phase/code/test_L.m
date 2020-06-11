function [probs1,probs2] = test_L(m)
n = 21;
Ls1 = linspace(2,5,n);
probs1 = [];
for i = 1:n
    L = Ls1(i);
    i
    count = 0;
    for j = 1:100
        err0 = phase_retrieval(m, 'gaussian', L, 0);
        if err0 < 1e-5
            count = count + 1;
        end
    end
    prob = count/100;
    probs1 = [probs1;prob];
end

n = 11;
Ls2 = linspace(2,12,n);
probs2 = [];
for i = 1:n
    L = Ls2(i);
    i
    count = 0;
    for j = 1:100
        err0 = phase_retrieval(m, 'cdp', L, 0);
        if err0 < 1e-5
            count = count + 1;
        end
    end
    prob = count/100;
    probs2 = [probs2; prob];
end


figure
hold on
subplot(1,2,1)
plot(Ls1, probs1,'-x', 'linewidth',1)
title('Gaussian model')
set(gca,'FontSize',12);
xlabel('$L$','interpreter','latex')
ylabel('probability of success')
grid on
subplot(1,2,2)
plot(Ls2, probs2,'-x', 'linewidth',1)
title('Coded diffusion model')
set(gca,'FontSize',12);
xlabel('$L$','interpreter','latex')
ylabel('probability of success')
grid on
