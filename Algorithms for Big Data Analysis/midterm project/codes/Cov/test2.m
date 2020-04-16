clear all; 
n = 30;
rhos = [0.001,0.01,0.1,1,10];
model = 2; % Choose one model to generate S

if model == 1
    S = model1(n);
else
    [S,flag] = model2(n);
    if flag ~= 1
        error('model 2 has no solution, please retry');
    end
end
S_inv = inv(S); 
results1 = zeros(length(rhos),6);
results2 = zeros(length(rhos),6);

for i = 1:length(rhos)
    rho = rhos(i);
    sigma = 1/rho
    tic;
    [X1, out1] = cvx_sol_2(S, sigma);
    out1
    t1 = toc;
    tic;
    [X2, out2] = spinv_admm(S, sigma);
    out2
    t2 = toc;
    results1(i,1) = sigma;
    results1(i,2) = out1.value;
    results1(i,3) = out1.dualgap;
    results1(i,4) = t1;
    results1(i,5) = nan;
    results2(i,1) = sigma;
    results2(i,2) = out2.value;
    results2(i,3) = out2.dualgap;
    results2(i,4) = t2;
    results2(i,5) = out2.iter;
    results2(i,6) = norm(X1-X2,'fro')/max(norm(X1,'fro'),1);
    name = 'extra_'+string(rho)+'_'+string(model)+'.eps';
    fig = figure;
    subplot(1,3,1)
    imagesc(abs(S_inv))
    colorbar
    map=colormap(flipud(gray));
    subplot(1,3,2)
    imagesc(abs(X1))
    colorbar
    subplot(1,3,3)
    imagesc(abs(X2))
    colorbar
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 1.6])
%     print(fig,name,'-depsc2')
end
