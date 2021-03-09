function adapt_test(iter, alpha, choice, err_type, derivative)
n0 = 50;
x = zeros(n0, 1);
for i = 1 : n0
    x(i) = i / n0;
end
if choice == 1
    g = 10*pi;
elseif choice == 2
    g = - 10 * exp(-2.5);
else
    error('Wrong choice');
end
N = zeros(iter, 1); ERR = zeros(iter, 1);
if choice == 1
    u = fem(@f1, g, x, 'adaptive');
else
    u = fem(@f2, g, x, 'adaptive');
end
for j = 1:iter
    if choice == 1
        post = post_err(@f1, x, u);
    else
        post = post_err(@f2, x, u);
    end
    threshold = alpha * max(post);
    for i = 1:length(u)
        if post(i) > threshold
            if i == 1
                x = [x' x(1)/2]';
            else
                x = [x' (x(i-1)+x(i))/2]';
            end
        end
    end
    x = sort(x);
    n = length(x); N(j) = n;
    if choice == 1
        u = fem(@f1, g, x, 'adaptive');
        real = sol1(x); dreal = dsol1(x);
    else
        u = fem(@f2, g, x, 'adaptive');
        real = sol2(x); dreal = dsol2(x);
    end
    if derivative
        h = zeros(n, 1);
        du = zeros(n, 1); du(1) = u(1)/x(1);
        for i = 2:n
            h(i) = x(i) - x(i-1);
            du(i) = (u(i) - u(i-1))/h(i);
        end
        Err = du - dreal;
        if strcmp(err_type, '2')
            err = sqrt(sum(Err.^2.*h));
            ERR(j) = err;
        elseif strcmp(err_type, 'inf')
            err = norm(Err, inf);
            ERR(j) = err;
        else
            error('wrong err type');
        end
        fig = figure;
        name = strcat('derivative\_case',string(choice),'\_', 'adaptive', '\_', string(length(x)),'.jpg');
        scatter(x, Err);
        hold on
        xlabel('x')
        ylabel('error')
        title(name)
        hold off
        img = frame2im(getframe(fig));
        name = strcat('derivative_case',string(choice),'_', 'adaptive', '_', string(length(x)),'.jpg');
        imwrite(img, name);
    else
        h = zeros(n, 1);
        for i = 2:n
            h(i) = x(i) - x(i-1);
        end
        Err = u - real;
        if strcmp(err_type, '2')
            err = sqrt(sum(Err.^2.*h));
            ERR(j) = err;
        elseif strcmp(err_type, 'inf')
            err = norm(Err, inf);
            ERR(j) = err;
        else
            error('wrong err type');
        end
        fig = figure;
        name = strcat('case',string(choice),'\_', 'adaptive', '\_', string(length(x)),'.jpg');
        scatter(x, Err);
        hold on
        xlabel('x')
        ylabel('error')
        title(name)
        hold off
        img = frame2im(getframe(fig));
        name = strcat('case',string(choice),'_', 'adaptive', '_', string(length(x)),'.jpg');
        imwrite(img, name);
    end
end
if derivative
    di = 'derivative\_';
else
    di = '';
end
log_err = log(ERR)
k = polyfit(log(N), log_err, 1);
if strcmp(err_type, 'inf')
    name = strcat(di, 'l^\infty err','\_','case',string(choice),'\_', 'adaptive','.jpg');
elseif strcmp(err_type, '2')
    name = strcat(di, 'l^2 err','\_','case',string(choice),'\_', 'adaptive','.jpg');
else
    err('wrong err type');
end
fig = figure;
loglog(N, ERR)
hold on
scatter(N,ERR)
xlabel('n');
ylabel('error');
title([name,'Slope: ',string(k(1))])
hold off
img = frame2im(getframe(fig));
if derivative
    di = 'derivative_';
else
    di = '';
end
if strcmp(err_type, 'inf')
    name = strcat(di, 'l infty err','_','case',string(choice),'_', 'adaptive','.jpg');
elseif strcmp(err_type, '2')
    name = strcat(di, 'l^2 err','_','case',string(choice),'_', 'adaptive','.jpg');
else
    err('wrong err type');
end
imwrite(img, name);
end

function post = post_err(f, x, u) 
n = length(x); post = zeros(n,1);
for i = 1:n
    if i == 1
        mean_f = (f(0) + f(x(1))) * x(1) / 2;
        post(1) = abs((x(1))^1.5 * sqrt(((u(1) - mean_f)^3 + (mean_f)^3) / 3 / u(1)));
    else
        mean_f = (f(x(i-1)) + f(x(i))) * (x(i) - x(i-1)) / 2;
        post(i) = abs((x(i) - x(i-1))^1.5 * sqrt(((u(i) - mean_f)^3 - (u(i-1) - mean_f)^3) / 3 / (u(i) - u(i-1))));
    end
end
end