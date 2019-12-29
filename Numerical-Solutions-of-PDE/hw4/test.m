% Calculate the error in one specific n. 
function [err,n] = test(n, choice, grid, err_type, derivative)
if strcmp(grid, 'uniform')
    x = zeros(n, 1);
    for i = 1 : n
        x(i) = i / n;
    end
elseif strcmp(grid, 'random')
    x = zeros(n, 1); x(1) = rand();
    for i = 2 : n
        x(i) = x(i-1) + rand();
    end
    x = x / x(n);
elseif strcmp(grid, 'perturbed')
    x = zeros(n, 1); x(n) = 1;
    rand('seed',0);
    for i = 1 : n-1
        x(i) = i/n  + (rand()-1/2) / 5 / n;
    end
else
    error('Wrong grid type');
end

n = length(x);
if choice == 1
    g = 10*pi; real = zeros(length(x), 1); dreal = zeros(length(x), 1);
    u = fem(@f1, g, x, grid);
    for i = 1 : length(x)
        real(i) = sol1(x(i));
        dreal(i) = dsol1(x(i));
    end
elseif choice == 2
    g = - 10 * exp(-2.5); real = zeros(length(x), 1); dreal = zeros(length(x), 1);
    u = fem(@f2, g, x, grid);
    for i = 1 : length(x)
        real(i) = sol2(x(i));
        dreal(i) = dsol2(x(i));
    end
else
    error('Wrong function choice');
end
if ~derivative
    h = zeros(n, 1); h(1) = x(1);
    for i = 2:n
        h(i) = x(i)-x(i-1);
    end
    Err = u - real;
    if strcmp(err_type, '2')
        err = sqrt(sum(Err.^2.*h));
    elseif strcmp(err_type, 'inf')
        err = norm(Err, inf);
    else
        error('wrong err type');
    end
    fig = figure;
    name = strcat('case',string(choice),'\_', grid, '\_', string(length(x)),'.jpg');
    scatter(x, Err);
    hold on
    xlabel('x')
    ylabel('error')
    title(name)
    hold off
    img = frame2im(getframe(fig));
    name = strcat('case',string(choice),'_', grid, '_', string(length(x)),'.jpg');
    imwrite(img, name);
else
    h = zeros(n, 1); h(1) = x(1);
    du = zeros(n, 1); du(1) = u(1)/h(1);
    for i = 2:n
        h(i) = x(i) - x(i-1);
        du(i) = (u(i) - u(i-1))/h(i);
    end
    Err = du - dreal;
    if strcmp(err_type, '2')
        err = sqrt(sum(Err.^2.*h));
    elseif strcmp(err_type, 'inf')
        err = norm(Err, inf);
    else
        error('wrong err type');
    end
    fig = figure;
    name = strcat('derivative\_case',string(choice),'\_', grid, '\_', string(length(x)),'.jpg');
    scatter(x, Err);
    hold on
    xlabel('x')
    ylabel('error')
    title(name)
    hold off
    img = frame2im(getframe(fig));
    name = strcat('derivative_case',string(choice),'_', grid, '_', string(length(x)),'.jpg');
    imwrite(img, name);
end
end