function main(choice)
if strcmp(choice, 'linear-discontinous-periodic')
    y1 = linear_periodic_solver(-1, 1, 500, 2500, 0.0032, 'upwind', @a, @init_value_discontinous_periodic);
    y2 = linear_periodic_solver(-1, 1, 500, 2500, 0.0032, 'LF', @a, @init_value_discontinous_periodic);
    y3 = linear_periodic_solver(-1, 1, 500, 2500, 0.0032, 'LW', @a, @init_value_discontinous_periodic);
    y4 = linear_periodic_solver(-1, 1, 500, 2500, 0.0032, 'BW', @a, @init_value_discontinous_periodic);
    x = -1:0.004:1; y_real = zeros(501,1);
    for i = 1:501
        y_real(i) = init_value_discontinous_periodic(x(i)-8);
    end
elseif strcmp(choice, 'linear-continous-periodic')
    t0 = 16;
    y1 = linear_periodic_solver(-1, 1, 500, t0/0.0032, 0.0032, 'upwind', @a, @init_value_continous_periodic);
    y2 = linear_periodic_solver(-1, 1, 500, t0/0.0032, 0.0032, 'LF', @a, @init_value_continous_periodic);
    y3 = linear_periodic_solver(-1, 1, 500, t0/0.0032, 0.0032, 'LW', @a, @init_value_continous_periodic);
    y4 = linear_periodic_solver(-1, 1, 500, t0/0.0032, 0.0032, 'BW', @a, @init_value_continous_periodic);
    x = -1:0.004:1; y_real = zeros(501,1);
    for i = 1:501
        y_real(i) = init_value_continous_periodic(x(i)-16);
    end
elseif strcmp(choice, 'linear-discontinous-infinite')
    n = 1/0.0025;
    y1 = linear_infinite_solver(0, 1, n, n, 1/2/n, 'upwind',@a, @init_value_discontinous_infinite);
    y2 = linear_infinite_solver(0, 1, n, n, 1/2/n, 'LF', @a, @init_value_discontinous_infinite);
    y3 = linear_infinite_solver(0, 1, n, n, 1/2/n, 'LW', @a, @init_value_discontinous_infinite);
    y4 = linear_infinite_solver(0, 1, n, n, 1/2/n, 'BW', @a, @init_value_discontinous_infinite);
    x = zeros(n+1, 1);
    for i = 1:n+1
        x(i) = (i-1)/n ;
    end
    y_real = zeros(n+1,1);
    for i = 1:n+1
        y_real(i) = init_value_discontinous_infinite(x(i)-0.5);
    end
elseif strcmp(choice, 'linear-continous-infinite')
    t = 8;
    y1 = linear_infinite_solver(t, t+1, 200, t/0.004, 0.004, 'upwind',@a, @init_value_continous_infinite);
    y2 = linear_infinite_solver(t, t+1, 200, t/0.004, 0.004, 'LF',@a, @init_value_continous_infinite);
    y3 = linear_infinite_solver(t, t+1, 200, t/0.004, 0.004, 'LW',@a, @init_value_continous_infinite);
    y4 = linear_infinite_solver(t, t+1, 200, t/0.004, 0.004, 'BW',@a, @init_value_continous_infinite);
    x = zeros(201, 1);
    for i = 1:201
        x(i) = (i-1)/200 + t ;
    end
    y_real = zeros(201,1);
    for i = 1:201
        y_real(i) = init_value_continous_infinite(x(i)-t);
    end
end
hold off
subplot(2,2,1)
hold on
scatter(x,y1,'red');
plot(x,y_real,'LineWidth',3,'Color','black');
title('Upwind')
hold off
subplot(2,2,2)
hold on
scatter(x,y2,'red');
plot(x,y_real,'LineWidth',3,'Color','black');
title('Lax-Friedrichs')
hold off
subplot(2,2,3)
hold on
scatter(x,y3,'red');
plot(x,y_real,'LineWidth',3,'Color','black');
title('Lax-Wendroff')
hold off
subplot(2,2,4)
title('Beam-Warming')
hold on
scatter(x,y4,'red');
plot(x,y_real,'LineWidth',3,'Color','black');
hold off
end
