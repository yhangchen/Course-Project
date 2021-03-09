%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This program implement numerical schemes of nonlinear advection equation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = nonlinear_solver(x_left, x_right, n, m, dt, scheme_choice, f, df, init_value)
x = zeros(n+2*m, 1); init = zeros(n+2*m, 1); h = (x_right - x_left) / n;
for i = 1 : 2 * m + n
    x(i) = x_left + (i - 0.5 - m) * h;
end
% Simpson integration
for i = 1 : 2 * m + n
    init(i) = (init_value(x(i) - h / 2) + init_value(x(i) + h / 2)...
        + 4 * init_value(x(i))) / 6;
end
if strcmp(scheme_choice, 'upwind')
    for j = 1 : m
        u = zeros(n+2*(m-j), 1);
        for i = 1 : n+2*(m-j)
            if abs(df(init(i+1)) * dt / h) > 1
                error('Disobeying CFL')
            end
            if init(i+1) == init(i+2)
                a1 = 0;
            else
                a1 = ((f(init(i+2))) - (f(init(i+1)))) / (init(i+2) - init(i+1));
            end
            if init(i+1) == init(i)
                a0 = 0;
            else
                a0 = ((f(init(i+1))) - (f(init(i)))) / (init(i+1) - init(i));
            end
            u(i) = init(i+1) - dt / 2 / h * (((1 + sign(a1)) * f(init(i+1)) ...
                + ((1 - sign(a1))) * f(init(i+2))) - ((1 + sign(a0)) * f(init(i))...
                + (1 - sign(a0)) * f(init(i+1))));
        end
        init = u;
    end
elseif strcmp(scheme_choice, 'Roe')
    for j = 1 : m
        u = zeros(n+2*(m-j), 1);
        for i = 1 : n+2*(m-j)
            if abs(df(init(i+1)) * dt / h) > 1
                error('Disobeying CFL')
            end
            if init(i+1) == init(i+2)
                a1 = 0;
            else
                a1 = ((f(init(i+2))) - (f(init(i+1)))) / (init(i+2) - init(i+1));
            end
            if init(i+1) == init(i)
                a0 = 0;
            else
                a0 = ((f(init(i+1))) - (f(init(i)))) / (init(i+1) - init(i));
            end
            u(i) = init(i+1) - dt / 2 / h * (f(init(i+2))-f(init(i)) - abs(a1)...
                * (init(i+2) - init(i+1)) + abs(a0) * (init(i+1) - init(i)));
        end
        init = u;
    end  
    
elseif strcmp(scheme_choice, 'Huang')
    for j = 1 : m
        u = zeros(n+2*(m-j), 1);
        for i = 1 : n+2*(m-j)
            if abs(df(init(i+1)) * dt / h) > 1
                error('Disobeying CFL')
            end
            u(i) = init(i+1) - dt / 2 / h * (f(init(i+2))-f(init(i)) - ...
                sign(df((init(i+2) + init(i+1))/2)) * (f(init(i+2)) - f(init(i+1)))...
                + sign(df((init(i+1) + init(i))/2)) * (f(init(i+1)) - f(init(i))));
        end
        init = u;
    end
elseif strcmp(scheme_choice, 'Osher')
    for j = 1 : m
        u = zeros(n+2*(m-j), 1);
        for i = 1 : n+2*(m-j)
            if abs(df(init(i+1)) * dt / h) > 1
                error('Disobeying CFL')
            end
            a1 = (init(i+2) - init(i+1)) * (abs(df(init(i+1))) + abs(df(init(i+2))) + 4 * abs(df((init(i+1)+init(i+2))/2)))/6;
            a0 = (init(i+1) - init(i)) * (abs(df(init(i+1))) + abs(df(init(i))) + 4 * abs(df((init(i+1)+init(i))/2)))/6;
            u(i) = init(i+1) - dt / 2 / h * (f(init(i+2))-f(init(i)) - a1 + a0);
        end
        init = u;
    end
    
elseif strcmp(scheme_choice, 'Richtmyer')
    for j = 1 : m
        u = zeros(n+2*(m-j), 1);
        for i = 1 : n+2*(m-j)
            if abs(df(init(i+1)) * dt / h) > 1
                error('Disobeying CFL')
            end
            u_middle = (init(i) + init(i+1)) / 2 - dt / 2 / h * (f(init(i+1)) - f(init(i)));
            u_new_middle = (init(i+1) + init(i+2)) / 2 - dt / 2 / h * (f(init(i+2)) - f(init(i+1)));
            u(i) = init(i+1) - dt / h * (f(u_new_middle) - f(u_middle));
        end
        init = u;
    end
elseif strcmp(scheme_choice, 'LW')
    for j = 1 : m
        u = zeros(n+2*(m-j), 1);
        for i = 1 : n+2*(m-j)
            if abs(df(init(i+1)) * dt / h) > 1
                error('Disobeying CFL')
            end
            u(i) = init(i+1) - dt / 2 / h * (f(init(i+2)) - f(init(i)))...
                + dt^2 / 2 / h^2 * (df((init(i+2) + init(i+1)) / 2) * (f(init(i+2)) - f(init(i+1))) ...
                - df((init(i+1) + init(i)) / 2) * (f(init(i+1)) - f(init(i))));
        end
        init = u;
    end
elseif strcmp(scheme_choice, 'MacCormack')
    for j = 1 : m
        u = zeros(n+2*(m-j), 1);
        for i = 1 : n+2*(m-j)
            if abs(df(init(i+1)) * dt / h) > 1
                error('Disobeying CFL')
            end
            u_star_2 = init(i+1) - dt / h * (f(init(i+2)) - f(init(i+1)));
            u_star_1 = init(i) - dt / h * (f(init(i+1)) - f(init(i)));
            u(i) = (init(i+1) + u_star_2) / 2 - dt / h / 2 * (f(u_star_2) - f(u_star_1));
        end
        init = u;
    end    
    
end
end