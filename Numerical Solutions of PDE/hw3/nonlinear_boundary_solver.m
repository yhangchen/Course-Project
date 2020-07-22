function u = nonlinear_boundary_solver(x_left, x_right, n, m, dt, scheme_choice, a, f, df, init_value, left, right)
x = zeros(n+2, 1); init = zeros(n+1, 1); h = (x_right - x_left) / n;
for i = 0 : n + 1
    x(i+1) = x_left + (i - 0.5) * h;
end
% Simpson integration
for i = 1 : n
    init(i) = (init_value(x(i+1) - h / 2) + init_value(x(i+1) + h / 2)...
        + 4 * init_value(x(i+1)))  / 6;
end
init(n+1) = 2*right(0) - init(n);

if strcmp(scheme_choice, 'upwind')
    for j = 1 : m
        u = zeros(n+1, 1);
        for i = 1 : n
            if abs(df(init(i)) * dt / h) > 1
                error('Disobeying CFL')
            end
            if init(i) == init(i+1)
                a1 = 0;
            else
                a1 = ((f(init(i+1))) - (f(init(i)))) / (init(i+1) - init(i));
            end
            if a1 >= 0
                F1 = f(init(i));
            else
                F1 = f(init(i+1));
            end
%             else
%                 a1 = df(right((j-1)*dt));
%                 if a1 >= 0
%                     F1 = (f(right(j*dt))+f(right((j-1)*dt))+4*f(right((j-1/2)*dt))) / 6;
%                 else
%                     F1 = f(init(n));
%                 end
%             end
            if i > 1
                if init(i) == init(i-1)
                    a0 = 0;
                else
                    a0 = (f(init(i)) - f(init(i-1))) / (init(i) - init(i-1));
                end
                if a0 >= 0
                    F0 = f(init(i-1));
                else
                    F0 = f(init(i));
                end
            else
                a0 = df(left((j-1)*dt));
                if a0 >= 0
                    F0 = (f(left(j*dt))+f(left((j-1)*dt))+4*f(left((j-1/2)*dt))) / 6;
                else
                    F0 = f(init(1));
                end
            end
            
            u(i) = init(i) - dt / h * (F1 - F0);
        end
        u(n+1) = 2*right(j*dt) - u(n);
        init = u;
    end
    u = u(1:n);
elseif strcmp(scheme_choice, 'Richtmyer middle')
    for j = 1 : m
        u = zeros(n, 1);
        for i = 1 : n
            if abs(df(init(i)) * dt / h) > 1
                error('Disobeying CFL')
            end
            if i > 1
                u_middle = (init(i-1) + init(i)) / 2 - dt / 2 / h * (f(init(i)) - f(init(i-1)));
            else
                u_middle = left((j-1/2) * dt);
            end
            if i < n
                u_new_middle = (init(i+1) + init(i)) / 2 - dt / 2 / h * (f(init(i+1)) - f(init(i)));
            else
                u_new_middle = right((j-1/2) * dt);
            end
            u(i) = init(i) - dt / h * (f(u_new_middle) - f(u_middle));
        end
        init = u;
    end
elseif strcmp(scheme_choice, 'Richtmyer Thompson')
    for j = 1 : m
        u = zeros(n, 1);
        for i = 1 : n
            if abs(df(init(i)) * dt / h) > 1
                error('Disobeying CFL')
            end
            if i > 1
                u_middle = f((init(i-1) + init(i)) / 2 - dt / 2 / h * (f(init(i)) - f(init(i-1))));
            else
                u_middle = (f(left(j*dt))+f(left((j-1)*dt))+4*f(left((j-1/2)*dt))) / 6;
            end
            if i < n
                u_new_middle = f((init(i+1) + init(i)) / 2 - dt / 2 / h * (f(init(i+1)) - f(init(i))));
            else
                u_new_middle = (f(right(j*dt))+f(right((j-1)*dt))+4*f(right((j-1/2)*dt))) / 6;
            end
            u(i) = init(i) - dt / h * (u_new_middle - u_middle);
        end
        init = u;
    end
else
    error('Wrong input')
end
end

