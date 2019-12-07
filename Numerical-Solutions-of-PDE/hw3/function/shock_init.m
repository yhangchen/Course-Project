function u0 = shock_init(x)
if x < 0
    u0 = 2;
else
    u0 = 1;
end
end