function err_nonlinear_finite(scheme_choice)
hold off
err = zeros(10,1); h = zeros(10,1); err2 = zeros(10,1); in_err = zeros(10,1);in_err2 = zeros(10,1);
subplot(6,2,[1,2])
for i = 1:10
    n = 100*i; m = 2*n; h(i) = 1/n;
    x = zeros(n,1); y_real = zeros(n,1);
    for j = 1:n
        x(j) = (j-1/2)/n;
        y_real(j) = real_nonlinear_boundary(x(j),1);
    end
    Err = nonlinear_boundary_solver(0, 1, n, m, 1/m, scheme_choice, @a, @f1, @df1, @init_value_nonlinear_finite, @left_nonlinear, @right_nonlinear) - y_real;
    plot(x,Err)
    hold on
    err(i) = norm(Err, inf);
    err2(i) = norm(Err);
    m = n/10;
    inner_Err = Err(m:n-m);
    in_err(i) = norm(inner_Err, inf);
    in_err2(i) = norm(inner_Err);
end
log_h = log(h);
xlabel('x');
ylabel('error');
title([scheme_choice,' t = 1']);
legend('n=100','n=200','n=300','n=400','n=500','n=600','n=700','n=800','n=900','n=1000')
hold off
log_err = log(err);
log_err2 = log(err2);
log_in_err = log(in_err);
log_in_err2 = log(in_err2);
k1 = polyfit(log_h, log_err, 1);
k2 = polyfit(log_h, log_err2, 1);
ik1 = polyfit(log_h, log_in_err, 1);
ik2 = polyfit(log_h, log_in_err2, 1);
subplot(6,2,3)
loglog(h, err)
hold on
scatter(h,err)
xlabel('h');
ylabel('error');
title(['l^\infty, Slope: ',string(k1(1))])
hold off
subplot(6,2,4)
loglog(h, err2)
hold on
scatter(h,err2)
xlabel('h');
ylabel('error');
title(['l^2, Slope: ',string(k2(1))])
hold off
subplot(6,2,5)
loglog(h, in_err)
hold on
scatter(h,in_err)
xlabel('h');
ylabel('error in [0.1,0.9]');
title(['l^\infty, Slope: ',string(ik1(1))])
hold off
subplot(6,2,6)
loglog(h, in_err2)
hold on
scatter(h,in_err2)
xlabel('h');
ylabel('error in [0.1,0.9]');
title(['l^2, Slope: ',string(ik2(1))])
hold off




subplot(6,2,[7,8])
for i = 1:10
    n = 100*i; m = 2*n; h(i) = 1/n;
    x = zeros(n,1); y_real = zeros(n,1);
    for j = 1:n
        x(j) = (j-1/2)/n;
        y_real(j) = real_nonlinear_boundary(x(j),1.5);
    end
    Err = nonlinear_boundary_solver(0, 1, n, m, 1.5/m, scheme_choice, @a, @f1, @df1, @init_value_nonlinear_finite, @left_nonlinear, @right_nonlinear) - y_real;
    plot(x,Err)
    hold on
    err(i) = norm(Err, inf);
    err2(i) = norm(Err);
end
xlabel('x');
ylabel('error');
title([scheme_choice, ' t = 1.5']);
legend('n=100','n=200','n=300','n=400','n=500','n=600','n=700','n=800','n=900','n=1000')
log_err = log(err);
log_err2 = log(err2);
log_in_err = log(in_err);
log_in_err2 = log(in_err2);
k1 = polyfit(log_h, log_err, 1);
k2 = polyfit(log_h, log_err2, 1);
ik1 = polyfit(log_h, log_in_err, 1);
ik2 = polyfit(log_h, log_in_err2, 1);
subplot(6,2,9)
loglog(h, err)
hold on
scatter(h,err)
xlabel('h');
ylabel('error');
title(['l^\infty, Slope: ',string(k1(1))])
hold off
subplot(6,2,10)
loglog(h, err2)
hold on
scatter(h,err2)
xlabel('h');
ylabel('error');
title(['l^2, Slope: ',string(k2(1))])
hold off
subplot(6,2,11)
loglog(h, in_err)
hold on
scatter(h,in_err)
xlabel('h');
ylabel('error in [0.1,0.9]');
title(['l^\infty, Slope: ',string(ik1(1))])
hold off
subplot(6,2,12)
loglog(h, in_err2)
hold on
scatter(h,in_err2)
xlabel('h');
ylabel('error in [0.1,0.9]');
title(['l^2, Slope: ',string(ik2(1))])
hold off
