function Err = convergence_test(choice)
if strcmp(choice, 'linear')
    Err1 = zeros(5, 5); Err2 = zeros(5, 5);
    Err3 = zeros(5, 5); Err4 = zeros(5, 5);
    for n0 = 1:5
        for t0 = 1:5
        n = 100 * n0; dt = 1.6 / n; t = 1.6 * t0;
        y1 = linear_periodic_solver(-1, 1, n, t/dt, dt, 'upwind', @a, @init_value_continous_periodic);
        y2 = linear_periodic_solver(-1, 1, n, t/dt, dt, 'LF', @a, @init_value_continous_periodic);
        y3 = linear_periodic_solver(-1, 1, n, t/dt, dt, 'LW', @a, @init_value_continous_periodic);
        y4 = linear_periodic_solver(-1, 1, n, t/dt, dt, 'BW', @a, @init_value_continous_periodic);
        x = -1:2/n:1; y_real = zeros(n+1,1);
        for i = 1:n+1
            y_real(i) = init_value_continous_periodic(x(i)-t);
        end
        Err1(n0,t0) = norm(y_real - y1, inf);
        Err2(n0,t0) = norm(y_real - y2, inf);
        Err3(n0,t0) = norm(y_real - y3, inf);
        Err4(n0,t0) = norm(y_real - y4, inf);
        end
    end
    Err = [Err1', Err2', Err3', Err4']';
elseif strcmp(choice, 'burgers')
    tlst = [0.8, 0.9];
    Err1 = zeros(5,2); Err2 = zeros(5,2);Err3 = zeros(5,2);Err4 = zeros(5,2);Err5 = zeros(5,2);Err6 = zeros(5,2);Err7 = zeros(5,2);
    for n0 = 1:5
        for t0 = 1:2
            t = tlst(t0); dt = 0.005/n0; n = 100*n0;
            y1 = nonlinear_solver(0, 2*pi, n, t/dt, dt, 'upwind',@f, @df, @init_value_continuous_burgers);
            y2 = nonlinear_solver(0, 2*pi, n, t/dt, dt, 'Richtmyer',@f, @df, @init_value_continuous_burgers);
            y3 = nonlinear_solver(0, 2*pi, n, t/dt, dt, 'LW',@f, @df, @init_value_continuous_burgers);
            y4 = nonlinear_solver(0, 2*pi, n, t/dt, dt, 'MacCormack',@f, @df, @init_value_continuous_burgers);
            y5 = nonlinear_solver(0, 2*pi, n, t/dt, dt, 'Roe',@f, @df, @init_value_continuous_burgers);
            y6 = nonlinear_solver(0, 2*pi, n, t/dt, dt, 'Huang',@f, @df, @init_value_continuous_burgers);
            y7 = nonlinear_solver(0, 2*pi, n, t/dt, dt, 'Osher',@f, @df, @init_value_continuous_burgers);
            x = zeros(n, 1); y_real = zeros(n, 1);
            for i = 1:n
                x(i) = (i-1/2)/n * 2 * pi;
                y_real(i) = newton_sol(x(i),t);
            end
            Err1(n0,t0) = norm(y_real - y1, inf);
            Err2(n0,t0) = norm(y_real - y2, inf);
            Err3(n0,t0) = norm(y_real - y3, inf);
            Err4(n0,t0) = norm(y_real - y4, inf);
            Err5(n0,t0) = norm(y_real - y5, inf);
            Err6(n0,t0) = norm(y_real - y6, inf);
            Err7(n0,t0) = norm(y_real - y7, inf);
        end
    end
    Err = [Err1', Err2', Err3', Err4', Err5', Err6', Err7'];
end
