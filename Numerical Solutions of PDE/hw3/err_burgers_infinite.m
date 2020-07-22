function err_burgers_infinite(scheme_choice)
hold off
err = zeros(10,1); h = zeros(10,1); err2 = zeros(10,1);
subplot(4,2,[1,2])
for i = 6:15
    n = 50*i+200; m = n; h(i-5) = 1/n;
    x = zeros(n,1); y_real = zeros(n,1);
    for j = 1:n
        x(j) = (j-1/2)/n * 2 * pi;
        y_real(j) = newton_sol(x(j),0.8);
    end
    Err = nonlinear_solver(0, 2*pi, n, m, 0.8/m, scheme_choice, @f, @df, @init_value_continuous_burgers) - y_real;
    plot(x,Err)
    hold on
    err(i-5) = norm(Err, inf);
    err2(i-5) = norm(Err);
end
log_h = log(h);
xlabel('x');
ylabel('error');
title([scheme_choice,' t = 0.8']);
legend('n=500','n=550','n=600','n=650','n=700','n=750','n=800','n=850','n=900','n=950')
hold off
log_err = log(err);
log_err2 = log(err2);
k1 = polyfit(log_h, log_err, 1);
k2 = polyfit(log_h, log_err2, 1);
subplot(4,2,3)
loglog(h, err)
hold on
scatter(h,err)
xlabel('h');
ylabel('error');
title(['Slope: ',string(k1(1))])
hold off
subplot(4,2,4)
loglog(h, err2)
hold on
scatter(h,err2)
xlabel('h');
ylabel('error');
title(['Slope: ',string(k2(1))])
hold off
subplot(4,2,[5,6])
for i = 6:15
    n = 50*i+200; m = 2*n; h(i-5) = 1/n;
    x = zeros(n,1); y_real = zeros(n,1);
    for j = 1:n
        x(j) = (j-1/2)/n * 2 * pi;
        y_real(j) = newton_sol(x(j),0.84);
    end
    Err = nonlinear_solver(0, 2*pi, n, m, 0.84/m, scheme_choice,@f, @df, @init_value_continuous_burgers) - y_real;
    plot(x,Err)
    hold on
    err(i-5) = norm(Err, inf);
    err2(i-5) = norm(Err);
end
xlabel('x');
ylabel('error');
title([scheme_choice, ' t = 0.84']);
legend('n=500','n=550','n=600','n=650','n=700','n=750','n=800','n=850','n=900','n=950')
hold off
log_err = log(err);
log_err2 = log(err2);
k1 = polyfit(log_h, log_err, 1);
k2 = polyfit(log_h, log_err2, 1);
subplot(4,2,7)
loglog(h, err)
hold on
scatter(h,err)
xlabel('h');
ylabel('error');
title(['Slope: ',string(k1(1))])
hold off
subplot(4,2,8)
loglog(h, err2)
hold on
scatter(h,err2)
xlabel('h');
ylabel('error');
title(['Slope: ',string(k2(1))])
hold off
