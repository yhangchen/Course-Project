function z = shock_real(x,t)
s = f(2) - f(1);
if x < s * t
    z = 2;
else
    z = 1;
end
end