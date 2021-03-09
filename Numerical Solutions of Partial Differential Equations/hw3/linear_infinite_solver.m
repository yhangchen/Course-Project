%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This program implement numerical schemes of linear advection equation on
% the real line.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = linear_infinite_solver(x_left, x_right, n, m, dt, scheme_choice, a, init_value)
x = zeros(2 * m + n + 1, 1); init = zeros(2 * m + n + 1, 1); h = (x_right - x_left) / n;
for i = 1 : 2 * m + n + 1
    x(i) = x_left + (i - m - 1) * h;
    init(i) = init_value(x(i));
end
% upwind scheme
if strcmp(scheme_choice, 'upwind')
    for j = 1:m
        u = zeros(n+1+2*(m-j), 1);
        for i = 1:n+1+2*(m-j)
            nu = a(x(i+j), j*dt) * dt / h;
            if abs(nu) > 1
                error('Disobeying CFL')
            end
            u(i) = init(i+1) - nu / 2 * (init(i+2) - init(i))...
                + abs(nu) / 2 * (init(i+2) + init(i) - 2 * init(i+1));           
        end
        init = u;
    end
% LF scheme
elseif strcmp(scheme_choice, 'LF')
    for j = 1:m    
        u = zeros(n+1+2*(m-j), 1);
        for i = 1:n+1+2*(m-j)
            nu = a(x(i+j), j*dt) * dt / h;
            if abs(nu) > 1
                error('Disobeying CFL')
            end
            u(i) = (init(i+2) + init(i)) / 2 - nu ...
                * (init(i+2) - init(i)) / 2;
        end
        init = u;
    end
elseif strcmp(scheme_choice, 'LW')
    for j = 1:m    
        u = zeros(n+1+2*(m-j), 1);
        for i = 1:n+1+2*(m-j)
            nu = a(x(i+j), j*dt) * dt / h;
            if abs(nu) > 1
                error('Disobeying CFL')
            end
            u(i) = -nu * (1- nu) / 2 * init(i+2) + (1 - nu^2) * init(i+1)...
                + nu * (1 + nu) / 2 * init(i);
        end
        init = u;
    end

% BW scheme
elseif strcmp(scheme_choice, 'BW')
    for j = 1:m
        u = zeros(n+1+2*(m-j), 1);
        for i = 1:n+1+2*(m-j)
            nu = a(x(i+j), j*dt) * dt / h;
            if abs(nu) > 2
                error('Disobeying CFL')
            end
            if nu > 0
                % on the boundary, we adopt LW scheme.
                if i == 1
                    u(i) = -nu * (1- nu) / 2 * init(i+2) + (1 - nu^2) * init(i+1)...
                        + nu * (1 + nu) / 2 * init(i);
                else
                    u(i) = (1 - nu) * (2 - nu) * init(i+1) / 2 + nu * (2 - nu) * init(i)...
                        - nu * (1 - nu) * init(i-1) / 2;
                end
            elseif nu <= 0
                if i == n+1+2*(m-j)
                    u(i) = -nu * (1- nu) / 2 * init(i+2) + (1 - nu^2) * init(i+1)...
                        + nu * (1 + nu) / 2 * init(i);
                else
                    u(i) = (1 + nu) * (2 + nu) * init(i+1) / 2 - nu * (2 + nu) * init(i+3)...
                        + nu * (1 + nu) * init(i+3) / 2;
                end
            end
        end
        init = u;
    end
end
end