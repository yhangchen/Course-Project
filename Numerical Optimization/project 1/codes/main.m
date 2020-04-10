function Result = main(choice)
% reproduce all results, choice = 1,2,3
if choice == 1
%% Brown Dennis function test
x0 = [25;5;-5;-1];
Result = [];
ms = [4,10,20,30,40,50];%4,10,20,30,40,50
for i = 1:length(ms)
    m = ms(i);
    fun = @(x) test_brown_den(x, m);

    opts = []; 
    opts.dnorm = 1e-2;
    tic
    [x1, out_vanilla] = newton_vanilla(fun, x0, opts);
    t = toc;
    out_vanilla
    results(1,1) = out_vanilla.value;
    results(2,1) = out_vanilla.gradnorm;
    results(3,1) = out_vanilla.iter;
    results(4,1) = out_vanilla.feva;
    results(5,1) = t;
    tic
    [x2, out_damped] = newton_damped(fun, x0, opts);
    t = toc;
    out_damped
    results(1,2) = out_damped.value;
    results(2,2) = out_damped.gradnorm;
    results(3,2) = out_damped.iter;
    results(4,2) = out_damped.feva;
    results(5,2) = t;
    
    tic
    [x3, out_mixed] = newton_mixed(fun, x0, opts);
    t = toc;
    out_mixed
    results(1,3) = out_mixed.value;
    results(2,3) = out_mixed.gradnorm;
    results(3,3) = out_mixed.iter;
    results(4,3) = out_mixed.feva;
    results(5,3) = t;
    
    tic
    [x4, out_LM] = newton_LM(fun, x0, opts);
    t = toc;
    out_LM
    results(1,4) = out_LM.value;
    results(2,4) = out_LM.gradnorm;
    results(3,4) = out_LM.iter;
    results(4,4) = out_LM.feva;
    results(5,4) = t;
    
    tic
    [x5, out_SR1] = newton_SR1(fun, x0, opts);
    t = toc;
    out_SR1
    results(1,5) = out_SR1.value;
    results(2,5) = out_SR1.gradnorm;
    results(3,5) = out_SR1.iter;
    results(4,5) = out_SR1.feva;
    results(5,5) = t;
    
    tic
    [x6, out_DFP] = newton_Broyden(fun, x0, 0, opts);
    t = toc;
    out_DFP
    results(1,6) = out_DFP.value;
    results(2,6) = out_DFP.gradnorm;
    results(3,6) = out_DFP.iter;
    results(4,6) = out_DFP.feva;
    results(5,6) = t;
    
    tic
    [x7, out_BFGS] = newton_Broyden(fun, x0, 1, opts);
    t = toc;
    out_BFGS
    results(1,7) = out_BFGS.value;
    results(2,7) = out_BFGS.gradnorm;
    results(3,7) = out_BFGS.iter;
    results(4,7) = out_BFGS.feva;
    results(5,7) = t;
    
    tic
    [x8, out_Broyden] = newton_Broyden(fun, x0, 0.5, opts);
    t = toc;
    out_Broyden
    results(1,8) = out_Broyden.value;
    results(2,8) = out_Broyden.gradnorm;
    results(3,8) = out_Broyden.iter;
    results(4,8) = out_Broyden.feva;
    results(5,8) = t;
    
    tic
    [x8, out_Broyden] = newton_Broyden(fun, x0, 2, opts);
    t = toc;
    out_Broyden
    results(1,9) = out_Broyden.value;
    results(2,9) = out_Broyden.gradnorm;
    results(3,9) = out_Broyden.iter;
    results(4,9) = out_Broyden.feva;
    results(5,9) = t;
    Result = [Result;results];



end
elseif choice == 2
%% discrete integral equation test
ns = [2, 10, 20, 30, 40, 50];
Result = [];
for i = 1:length(ns)
    n = ns(i);
    t = [1:n]'/(n+1);
    x0 = t.*(t-1);
    fun = @(x) test_disc_ie(x);
    opts = []; opts.dnorm = 1e-2;
    tic
    [x1, out_vanilla] = newton_vanilla(fun, x0, opts);
    t = toc;
    out_vanilla
    results(1,1) = out_vanilla.value;
    results(2,1) = out_vanilla.gradnorm;
    results(3,1) = out_vanilla.iter;
    results(4,1) = out_vanilla.feva;
    results(5,1) = t;
    tic
    [x2, out_damped] = newton_damped(fun, x0, opts);
    t = toc;
    out_damped
    results(1,2) = out_damped.value;
    results(2,2) = out_damped.gradnorm;
    results(3,2) = out_damped.iter;
    results(4,2) = out_damped.feva;
    results(5,2) = t;
    
    tic
    [x3, out_mixed] = newton_mixed(fun, x0, opts);
    t = toc;
    out_mixed
    results(1,3) = out_mixed.value;
    results(2,3) = out_mixed.gradnorm;
    results(3,3) = out_mixed.iter;
    results(4,3) = out_mixed.feva;
    results(5,3) = t;
    
    tic
    [x4, out_LM] = newton_LM(fun, x0, opts);
    t = toc;
    out_LM
    results(1,4) = out_LM.value;
    results(2,4) = out_LM.gradnorm;
    results(3,4) = out_LM.iter;
    results(4,4) = out_LM.feva;
    results(5,4) = t;
    
    tic
    [x5, out_SR1] = newton_SR1(fun, x0, opts);
    t = toc;
    out_SR1
    results(1,5) = out_SR1.value;
    results(2,5) = out_SR1.gradnorm;
    results(3,5) = out_SR1.iter;
    results(4,5) = out_SR1.feva;
    results(5,5) = t;
    
    tic
    [x6, out_DFP] = newton_Broyden(fun, x0, 0, opts);
    t = toc;
    out_DFP
    results(1,6) = out_DFP.value;
    results(2,6) = out_DFP.gradnorm;
    results(3,6) = out_DFP.iter;
    results(4,6) = out_DFP.feva;
    results(5,6) = t;
    
    tic
    [x7, out_BFGS] = newton_Broyden(fun, x0, 1, opts);
    t = toc;
    out_BFGS
    results(1,7) = out_BFGS.value;
    results(2,7) = out_BFGS.gradnorm;
    results(3,7) = out_BFGS.iter;
    results(4,7) = out_BFGS.feva;
    results(5,7) = t;
    
    tic
    [x8, out_Broyden] = newton_Broyden(fun, x0, 0.5, opts);
    t = toc;
    out_Broyden
    results(1,8) = out_Broyden.value;
    results(2,8) = out_Broyden.gradnorm;
    results(3,8) = out_Broyden.iter;
    results(4,8) = out_Broyden.feva;
    results(5,8) = t;
    
    tic
    [x8, out_Broyden] = newton_Broyden(fun, x0, 2, opts);
    t = toc;
    out_Broyden
    results(1,9) = out_Broyden.value;
    results(2,9) = out_Broyden.gradnorm;
    results(3,9) = out_Broyden.iter;
    results(4,9) = out_Broyden.feva;
    results(5,9) = t;
    Result = [Result;results];

end
elseif choice == 3
ns = [5,10,15,20,30,40];
Result = [];
for i = 1:length(ns)
    n = ns(i)
    hx = 1/n; hy = 1/n;
    xs = hx*[0:n] - 0.5;
    ys = hy*[0:n] - 0.5;
    x0 = zeros((n-1)^2,1);
    left = zeros(n+1,1); right = zeros(n+1,1);
    bottom = zeros(n+1,1); top = zeros(n+1,1);
    for j = 1:n+1
        left(j) = Enneper(-0.5,ys(j));
        right(j) = Enneper(0.5,ys(j));
        bottom(j) = Enneper(xs(j),-0.5);
        top(j) = Enneper(xs(j),0.5);
    end
    fun = @(x) test_minima_surface(x, n, n, bottom, top, left, right);
    
    opts = [];
    if n == 5
        opts.dnorm = 1e-3;
    else
        opts.dnorm = 1e-6;
    end
    
    tic
    [x5, out_SR1] = newton_SR1(fun, x0, opts);
    t = toc;
    out_SR1
    results(1,1) = out_SR1.value;
    results(2,1) = out_SR1.gradnorm;
    results(3,1) = out_SR1.iter;
    results(4,1) = out_SR1.feva;
    results(5,1) = t;
    
    tic
    [x6, out_DFP] = newton_Broyden(fun, x0, 0, opts);
    t = toc;
    out_DFP
    results(1,2) = out_DFP.value;
    results(2,2) = out_DFP.gradnorm;
    results(3,2) = out_DFP.iter;
    results(4,2) = out_DFP.feva;
    results(5,2) = t;
    
    tic
    [x7, out_BFGS] = newton_Broyden(fun, x0, 1, opts);
    t = toc;
    out_BFGS
    results(1,3) = out_BFGS.value;
    results(2,3) = out_BFGS.gradnorm;
    results(3,3) = out_BFGS.iter;
    results(4,3) = out_BFGS.feva;
    results(5,3) = t;
    
    tic
    [x8, out_Broyden] = newton_Broyden(fun, x0, 0.5, opts);
    t = toc;
    out_Broyden
    results(1,4) = out_Broyden.value;
    results(2,4) = out_Broyden.gradnorm;
    results(3,4) = out_Broyden.iter;
    results(4,4) = out_Broyden.feva;
    results(5,4) = t;
    
    tic
    [x9, out_Broyden_2] = newton_Broyden(fun, x0, 2, opts);
    t = toc;
    out_Broyden_2
    results(1,5) = out_Broyden_2.value;
    results(2,5) = out_Broyden_2.gradnorm;
    results(3,5) = out_Broyden_2.iter;
    results(4,5) = out_Broyden_2.feva;
    results(5,5) = t;
    Result = [Result;results];
    
    data = [];
    data.in = reshape(x7,n-1,n-1);
    data.left = left;
    data.right = right;
    data.bottom = bottom;
    data.top = top;
    Enneper_surface(n,n,data);
end

end
end