function u = newton_sol(x,t)
u = 0;
while abs(sin(x-u*t)+0.5-u) > 1e-6
    u = u - (sin(x-u*t)+0.5-u)/(cos(x-u*t)*(-t)-1);
end
