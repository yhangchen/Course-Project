% main function, output general results
function convergence(choice, grid, err_type, derivative)
if strcmp(grid, 'uniform') | strcmp(grid, 'random') | strcmp(grid, 'perturbed')
    if choice == 1
        s = 40:50;
    else
        s = 20:30;
    end
    n = 200*s;
    n = n';
else
    error('Wrong grid type');
end

m = length(n);
err = zeros(m, 1); n0 = zeros(m, 1);
for i = 1:m
    [err(i),n0(i)] = test(n(i), choice, grid, err_type, derivative);
end
log_err = log(err)
if derivative
    di = 'derivative\_';
else
    di = '';
end
k = polyfit(log(n), log_err, 1);
if strcmp(err_type, 'inf')
    name = strcat(di, 'l^\infty err','\_','case',string(choice),'\_', grid,'.jpg');
elseif strcmp(err_type, '2')
    name = strcat(di, 'l^2 err','\_','case',string(choice),'\_', grid,'.jpg');
else
    err('wrong err type');
end
fig = figure;
loglog(n0, err)
hold on
scatter(n0,err)
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
    name = strcat(di, 'l infty err','_','case',string(choice),'_', grid,'.jpg');
elseif strcmp(err_type, '2')
    name = strcat(di, 'l^2 err','_','case',string(choice),'_', grid,'.jpg');
else
    err('wrong err type');
end
imwrite(img, name);
end