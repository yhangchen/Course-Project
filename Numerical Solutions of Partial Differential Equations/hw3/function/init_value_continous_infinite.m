function u0 = init_value_continous_infinite(x,k)
k = 80;
u0 = exp(-100*(x-0.5)^2)*sin(k*x);
end