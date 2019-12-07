function u0 = init_value_discontinous_periodic(x)
x = mod(x+1, 2) - 1;
if -0.7 <= x <= 1
    xi = x - 0.3;
else
    xi = x + 1.7;
end
if -1 <= xi < -1/3
    u0 = -xi * sin(1.5*pi*xi^2);
elseif abs(xi) <= 1/3
    u0 = abs(sin(2*pi*xi));
else
    u0 = 2*xi - 1 - sin(3*pi*xi)/6;
end
end
