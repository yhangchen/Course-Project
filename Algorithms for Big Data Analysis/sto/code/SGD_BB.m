function [Value, Grad, Loss, w1, t] = SGD_BB(f, x, y, w, lambda, m, beta, lr0, lr1, max_iter, batchsize, fullsize, stopping)
[value1, grad1, loss1] = f(x, y, w, lambda);
Value = [value1]; Grad = [norm(grad1,inf)]; Loss = [loss1]; t = 0; 
c = 1; w1 = w; value1 = inf;
while t < max_iter
    if t == 0
        lr = lr0;
    elseif t == 1
        lr = lr1;
    else
        lr = norm(w1-w0)^2 / abs(dot(w1-w0, g1-g0)) / m;
        if isnan(lr) || isinf(lr)
            lr = lr1;
        end
        c = c^((t-2)/(t-1))*(lr*phi(t))^(1/(t-1));
        lr = c/phi(t);
    end
    
    if t > 0
        g0 = g1;
        w0 = w1;
    else
        g1 = 0;
    end
    loss11 = 0; value11 = 0; grad11 = 0;
    for i = 1:m
        ind = randi(fullsize, [1,batchsize]);
        x0 = x(:,ind);
        y0 = y(:,ind);
        [value10, grad, loss10] = f(x0, y0, w1, lambda);
        loss11 = loss11 + loss10; value11 = value11 + value10;
        grad11 = grad11 + norm(grad,inf);
        w1 = w1 - lr * grad;
        w1 = sign(w1) .* max(abs(w1) - lambda*lr, 0); % soft thresholding
        g1 = (1-beta) * g1 + beta * grad;
    end
    value0 = value1;
    value1 = value11/m;
    Loss = [Loss loss11/m];
    Grad = [Grad grad11/m];
    Value = [Value value11/m];
    t = t + 1;
    if abs(value0-value1) < stopping
        break
    end
end
end
function x = phi(t)
x = t;
end
    