function [value, grad, loss, w, m, v] = ADAM(f, x, y, w, lambda, m, v, t, lr, epsilon, beta1, beta2)
[~,grad,~] = f(x, y, w, lambda);
m = beta1 * m + (1 - beta1) * grad;
v = beta2 * v + (1 - beta2) * grad.^2;
m1 = m/(1 - beta1^t);
v1 = v/(1 - beta2^t);
w = w - lr * m1 ./ (sqrt(v1) + epsilon);
w = sign(w) .* max(abs(w) - lambda*lr ./ (sqrt(v1) + epsilon), 0);
[value, grad ,loss] = f(x, y, w, lambda);
end
