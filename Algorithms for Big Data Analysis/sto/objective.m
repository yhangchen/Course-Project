function [value, grad, bi_loss] = objective(x, y, w, lambda)
z = w'*x;
v = exp(-y.*z);
grad = (v.*(-y.*x))./(1+v);
value = mean(log(1+v))+ lambda*norm(w, 1);
bi_loss = mean((y.*z < 0));
grad = mean(grad, 2) + lambda*sign(w);
%mean(log(1+v)) + lambda*norm(w, 1);
end
