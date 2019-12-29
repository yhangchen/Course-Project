% derivative of u in case 2
function z = dsol2(x)
z = (10 - 20 * x) .* exp(-10*(x-0.5).^2);
end