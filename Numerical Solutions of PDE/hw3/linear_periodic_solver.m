%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This program implement numerical schemes under periodic boundary condition
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = linear_periodic_solver(x_left, x_right, n, m, dt, scheme_choice, a, initial)
x = zeros(n+1, 1); init = zeros(n+1, 1);
for i = 0:n
    x(i+1) = x_left + i * (x_right - x_left) / n;
    init(i+1) = initial(x(i+1));
end
if strcmp(scheme_choice, 'upwind')
    for j = 1:m
        u = zeros(n+1, 1);
        init(n+2) = init(2);
        for i = 2:n+1
            nu = a(x(i), j*dt) * dt * n / (x_right - x_left);
            if abs(nu) > 1
                error('Disobeying CFL')
            end
            u(i) = init(i) - nu / 2 * (init(i+1) - init(i-1))...
                + abs(nu) / 2 * (init(i+1) + init(i-1) - 2 * init(i));           
        end
        u(1) = u(n+1);
        init = u;
    end
elseif strcmp(scheme_choice, 'LF')
    for j = 1:m    
        u = zeros(n+1, 1);
        init(n+2) = init(2);
        for i = 2:n+1
            nu = a(x(i), j*dt) * dt * n / (x_right - x_left);
            if abs(nu) > 1
                error('Disobeying CFL')
            end
            u(i) = (init(i+1) + init(i-1)) / 2 - nu ...
                * (init(i+1) - init(i-1)) / 2;
        end
        u(1) = u(n+1);
        init = u;
    end
elseif strcmp(scheme_choice, 'LW')
    for j = 1:m    
        u = zeros(n+1, 1);
        init(n+2) = init(2);
        for i = 2:n+1
            nu = a(x(i), j*dt) * dt * n / (x_right - x_left);
            if abs(nu) > 1
                error('Disobeying CFL')
            end
            u(i) = -nu * (1- nu) / 2 * init(i+1) + (1 - nu^2) * init(i)...
                + nu * (1 + nu) / 2 * init(i-1);
        end
        u(1) = u(n+1);
        init = u;
    end
elseif strcmp(scheme_choice, 'BW')
    for j = 1:m
        u = zeros(n+1, 1);
        init(n+2 : n+4) = init(2 : 4); x(n+2) = x(2);
        for i = 3:n+2
            nu = a(x(i), j*dt) * dt * n / (x_right - x_left);
            if abs(nu) > 2
                error('Disobeying CFL')
            end
            if nu > 0
                u(i) = (1 - nu) * (2 - nu) * init(i) / 2 + nu * (2 - nu) * init(i-1)...
                    - nu * (1 - nu) * init(i-2) / 2;
            elseif nu <= 0
                u(i) = (1 + nu) * (2 + nu) * init(i) / 2 - nu * (2 + nu) * init(i+1)...
                    + nu * (1 + nu) * init(i+2) / 2;
            end
        end
        u(1) = u(n+1); u(2) = u(n+2);
        u = u(1: n+1);
        init = u;
    end
end
end
