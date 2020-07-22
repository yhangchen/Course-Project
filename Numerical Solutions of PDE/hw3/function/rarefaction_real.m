function z = rarefaction_real(x,t)
if x < t
    z = 1;
elseif (t <= x) & (x <= 2*t)
    z = x/t;
else
    z = 2;
end
end