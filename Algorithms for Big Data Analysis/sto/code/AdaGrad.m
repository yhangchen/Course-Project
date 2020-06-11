function [value, grad, loss, w, r] = AdaGrad(f, x, y, w, lambda, r, t, lr, epsilon)
[~, grad, ~] = f(x, y, w, lambda);
r = r + grad.^2;
w = w - lr *grad./ (sqrt(r) + epsilon);
w = sign(w) .* max(abs(w) - lambda*lr./ (sqrt(r) + epsilon), 0);
[value, grad, loss] = f(x, y, w, lambda);
end
