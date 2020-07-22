function main_burgers(choice)
if strcmp(choice,'continuous')
    t = 0.9; dt = 0.004; n = 200;
    y1 = nonlinear_solver(0, 2*pi, n, t/dt, dt, 'upwind',@f, @df, @init_value_continuous_burgers);
    y2 = nonlinear_solver(0, 2*pi, n, t/dt, dt, 'Richtmyer',@f, @df, @init_value_continuous_burgers);
    y3 = nonlinear_solver(0, 2*pi, n, t/dt, dt, 'LW',@f, @df, @init_value_continuous_burgers);
    y4 = nonlinear_solver(0, 2*pi, n, t/dt, dt, 'MacCormack',@f, @df, @init_value_continuous_burgers);
    y5 = nonlinear_solver(0, 2*pi, n, t/dt, dt, 'Roe',@f, @df, @init_value_continuous_burgers);
    y6 = nonlinear_solver(0, 2*pi, n, t/dt, dt, 'Huang',@f, @df, @init_value_continuous_burgers);
    y7 = nonlinear_solver(0, 2*pi, n, t/dt, dt, 'Osher',@f, @df, @init_value_continuous_burgers);
    x = zeros(n, 1);
    for i = 1:n
        x(i) = (i-1/2)/n * 2 * pi;
    end
    hold off
    subplot(1,3,1);
    hold on
    scatter(x,y1);
    scatter(x,y5);
    scatter(x,y6);
    scatter(x,y7);
    ylim([-1,2]);
    legend('Upwind', 'Roe','Huang','Engquist and Osher')
    title('Upwind type schemes');
    hold off
    subplot(1,3,2);
    hold on
    scatter(x,y2);
    scatter(x,y3);
    scatter(x,y4);
    ylim([-1,2]);
    legend('Richtmyer', 'Lax-Wendroff','MacCormack');
    title('Second-order schemes')
    hold off
    subplot(1,3,3);
    hold on
    scatter(x,y1);
    scatter(x,y5);
    scatter(x,y6);
    scatter(x,y7);
    scatter(x,y2);
    scatter(x,y3);
    scatter(x,y4);
    ylim([-1,2]);
    legend('Upwind','Roe','Huang','Engquist and Osher','Richtmyer', 'Lax-Wendroff','MacCormack')
    title('All the schemes')
    hold off

elseif strcmp(choice,'Riemann')
    t = 0.5; dt = 0.004; n = 400;
    y1 = nonlinear_solver(-2, 2, n, t/dt, dt, 'upwind',@f, @df, @rarefaction_init);
    y2 = nonlinear_solver(-2, 2, n, t/dt, dt, 'Richtmyer',@f, @df, @rarefaction_init);
    y3 = nonlinear_solver(-2, 2, n, t/dt, dt, 'LW',@f, @df, @rarefaction_init);
    y4 = nonlinear_solver(-2, 2, n, t/dt, dt, 'MacCormack',@f, @df, @rarefaction_init);
    y5 = nonlinear_solver(-2, 2, n, t/dt, dt, 'Roe',@f, @df, @rarefaction_init);
    y6 = nonlinear_solver(-2, 2, n, t/dt, dt, 'Huang',@f, @df, @rarefaction_init);
    y7 = nonlinear_solver(-2, 2, n, t/dt, dt, 'Osher',@f, @df, @rarefaction_init);
    
    z1 = nonlinear_solver(-2, 2, n, t/dt, dt, 'upwind',@f, @df, @shock_init);
    z2 = nonlinear_solver(-2, 2, n, t/dt, dt, 'Richtmyer',@f, @df, @shock_init);
    z3 = nonlinear_solver(-2, 2, n, t/dt, dt, 'LW',@f, @df, @shock_init);
    z4 = nonlinear_solver(-2, 2, n, t/dt, dt, 'MacCormack',@f, @df, @shock_init);
    z5 = nonlinear_solver(-2, 2, n, t/dt, dt, 'Roe',@f, @df, @shock_init);
    z6 = nonlinear_solver(-2, 2, n, t/dt, dt, 'Huang',@f, @df, @shock_init);
    z7 = nonlinear_solver(-2, 2, n, t/dt, dt, 'Osher',@f, @df, @shock_init);

    x = zeros(n, 1); y_real = zeros(n, 1); z_real = zeros(n, 1);
    for i = 1:n
        x(i) = (i-1/2)/n * 4 - 2;
        y_real(i) = rarefaction_real(x(i), t);
        z_real(i) = shock_real(x(i), t);
    end

    hold off
    subplot(2,3,1);
    hold on
    scatter(x,y1);
    scatter(x,y5);
    scatter(x,y6);
    scatter(x,y7);
    ylim([0.5,2.5]);
    plot(x,y_real,'LineWidth',1.5,'Color','black');
    legend('Upwind', 'Roe','Huang','Engquist and Osher')
    title('Upwind type schemes');
    hold off
    subplot(2,3,2);
    hold on
    scatter(x,y2);
    scatter(x,y3);
    scatter(x,y4);
    ylim([0.5,2.5]);
    plot(x,y_real,'LineWidth',1.5,'Color','black');
    legend('Richtmyer', 'Lax-Wendroff','MacCormack');
    title('Second-order schemes')
    hold off
    subplot(2,3,3);
    hold on
    scatter(x,y1);
    scatter(x,y5);
    scatter(x,y6);
    scatter(x,y7);
    scatter(x,y2);
    scatter(x,y3);
    scatter(x,y4);
    ylim([0.5,2.5]);
    plot(x,y_real,'LineWidth',1.5,'Color','black');
    legend('Upwind','Roe','Huang','Engquist and Osher','Richtmyer', 'Lax-Wendroff','MacCormack')
    title('All the schemes')
    hold off

    hold off
    subplot(2,3,4);
    hold on
    scatter(x,z1);
    scatter(x,z5);
    scatter(x,z6);
    scatter(x,z7);
    ylim([0.5,2.5]);
    plot(x,z_real,'LineWidth',1.5,'Color','black');
    legend('Upwind', 'Roe','Huang','Engquist and Osher')
    title('Upwind type schemes');
    hold off
    subplot(2,3,5);
    hold on
    scatter(x,z2);
    scatter(x,z3);
    scatter(x,z4);
    ylim([0.5,2.5]);
    plot(x,z_real,'LineWidth',1.5,'Color','black');
    legend('Richtmyer', 'Lax-Wendroff','MacCormack');
    title('Second-order schemes')
    hold off
    subplot(2,3,6);
    hold on
    scatter(x,z1);
    scatter(x,z5);
    scatter(x,z6);
    scatter(x,z7);
    scatter(x,z2);
    scatter(x,z3);
    scatter(x,z4);
    ylim([0.5,2.5]);
    plot(x,z_real,'LineWidth',1.5,'Color','black');
    legend('Upwind','Roe','Huang','Engquist and Osher','Richtmyer', 'Lax-Wendroff','MacCormack')
    title('All the schemes')
    hold off



end