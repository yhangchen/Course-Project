function u0 = rarefaction_init(x)
if x < 0
    u0 = 1;
else
    u0 = 2;
end
end