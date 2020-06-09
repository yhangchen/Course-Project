function [z, errs] = Wirtinger(A, At, y, x, max_iter, tau0)
% Initialization
m = length(y); n = length(x);
power_iter = 50; % Number of power iterations
z = randn(n,1); z = z/norm(z);
for i = 1:power_iter
    z = At(y.*A(z));
    z = z/norm(z);
end

lambda = sqrt(sum(y(:))/numel(y));
z = lambda * z;
normx = norm(x);
errs = norm(x - exp(-1i*angle(trace(x'*z))) * z)/normx;% Gradient Descent
mu = @(tau) min(1-exp(-tau/tau0), 0.2); % Schedule for step size
% gradient descent

for t = 1:max_iter
    w = A(z);
    grad  = At((abs(w).^2-y).*w); % Wirtinger gradient
    z = z - mu(t)/lambda^2 * grad;             % Gradient update 
    err = norm(x - exp(-1i*angle(x'*z)) * z)/normx;
%     if err < 1e-15
%         break
%     end
    errs = [errs, err];  
end
end

    
  
