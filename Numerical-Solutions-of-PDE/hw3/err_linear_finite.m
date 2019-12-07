function err_linear_finite(scheme_choice)
hold off
err = zeros(10,1); h = zeros(10,1); err2 = zeros(10,1);
subplot(2,2,[1,2])
for i = 1:10
    n = 100*i; m = n; h(i) = 1/n;
    x = linspace(0,1,n+1);
    real = zeros(n+1,1);
    for j = 1:n+1
        real(j) = real_linear_finite(x(j),1);
    end
    Err = linear_finite_solver(0, 1, n, m, 1/m, scheme_choice) - real;
    plot(x,Err)
    hold on
    err(i) = norm(Err, inf);
    err2(i) = norm(Err);
end
xlabel('x');
ylabel('error');
title(scheme_choice);
legend('n=100','n=200','n=300','n=400','n=500','n=600','n=700','n=800','n=900','n=1000')
subplot(2,2,3)
loglog(h, err)
hold on
scatter(h,err)
xlabel('h');
ylabel('error');
hold off
subplot(2,2,4)
loglog(h, err2)
hold on
scatter(h,err2)
xlabel('h');
ylabel('error');
hold off