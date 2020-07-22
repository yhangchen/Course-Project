%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This program implement numerical schemes on finite interval.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = linear_finite_solver(x_left, x_right, n, m, dt, scheme_choice)
x = zeros(n+1, 1); init = zeros(n+1, 1);
for i = 0:n
    x(i+1) = x_left + i * (x_right - x_left) / n;
    init(i+1) = init_value_linear_finite(x(i+1));
end

if strcmp(scheme_choice, 'upwind')
    for j = 1:m
        u = zeros(n+1, 1);
        u(1) = left_bd_linear_finite(j*dt); u(n+1) = right_bd_linear_finite(j*dt);
        for i = 2:n
            nu = a(x(i), j*dt) * dt * n / (x_right - x_left);
            if abs(nu) > 1
                error('Disobeying CFL')
            end
            u(i) = init(i) - nu / 2 * (init(i+1) - init(i-1))...
                + abs(nu) / 2 * (init(i+1) + init(i-1) - 2 * init(i));           
        end
        init = u;
    end
elseif strcmp(scheme_choice, 'LF')
    for j = 1:m    
        u = zeros(n+1, 1);
        u(1) = left_bd_linear_finite(j*dt); u(n+1) = right_bd_linear_finite(j*dt);
        for i = 2:n
            nu = a(x(i), j*dt) * dt * n / (x_right - x_left);
            if abs(nu) > 1
                error('Disobeying CFL')
            end
            u(i) = (init(i+1) + init(i-1)) / 2 - a(x(i), j*dt) * dt * n ...
                * (init(i+1) - init(i-1)) / 2;
        end
        init = u;
    end
elseif strcmp(scheme_choice, 'LW')
    for j = 1:m    
        u = zeros(n+1, 1);
        u(1) = left_bd_linear_finite(j*dt); u(n+1) = right_bd_linear_finite(j*dt);
        for i = 2:n
            nu = a(x(i), j*dt) * dt * n / (x_right - x_left);
            if abs(nu) > 1
                error('Disobeying CFL')
            end
            u(i) = -nu * (1- nu) / 2 * init(i+1) + (1 - nu^2) * init(i)...
                + nu * (1 + nu) / 2 * init(i-1);
        end
        init = u;
    end
elseif strcmp(scheme_choice, 'BW')
    for j = 1:m
        u = zeros(n+1, 1);
        u(1) = left_bd_linear_finite(j*dt); u(n+1) = right_bd_linear_finite(j*dt);
        for i = 2:n
            nu = a(x(i), j*dt) * dt * n / (x_right - x_left);
            if abs(nu) > 2
                error('Disobeying CFL')
            end
            if nu > 0
                % on the boundary, we adopt LW scheme.
                if i == 2
                    u(i) = -nu * (1- nu) / 2 * init(i+1) + (1 - nu^2) * init(i)...
                        + nu * (1 + nu) / 2 * init(i-1);
                else
                    u(i) = (1 - nu) * (2 - nu) * init(i) / 2 + nu * (2 - nu) * init(i-1)...
                        - nu * (1 - nu) * init(i-2) / 2;
                end
            elseif nu <= 0
                if i == n
                    u(i) = -nu * (1- nu) / 2 * init(i+1) + (1 - nu^2) * init(i)...
                        + nu * (1 + nu) / 2 * init(i-1);
                else
                    u(i) = (1 + nu) * (2 + nu) * init(i) / 2 - nu * (2 + nu) * init(i+1)...
                        + nu * (1 + nu) * init(i+2) / 2;
                end
            end
        end
        init = u;
    end
end
end