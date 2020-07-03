function lagrange(problem)
[nStations, nodes, nArcs, arc, nCommodities, commodities, cost]=readdata(problem);

B = zeros(nStations,nArcs);
for i = 1:nArcs
    B(arc(i,1),i) = 1;
    B(arc(i,2),i) = -1;
end
b = cell(1,nCommodities);
for k = 1:nCommodities
    btmp = zeros(nStations,1);
    for i = 1:nStations
        if commodities(k,1) == i
            btmp(i) = commodities(k,3);
        elseif commodities(k,2) == i
            btmp(i) = -commodities(k,3);
        else
            btmp(i) = 0;
        end
    b{k} = btmp;
    end
end
cap = arc(:,4);

UB = 1; LB = 0; tol = 1e-6;

iter = 0; flag = 0; lambda = 1; mu = ones(nArcs,1);

solution = cell(1,nCommodities);


maxIter = 100;
while iter < maxIter
    obj = 0;
    for k = 1:nCommodities
    [solution{k},objfunction,~,~,~]=linprog(mu,[],[],B,b{k},zeros(nArcs,1),[]);
    obj = obj + objfunction;
    end
    
    Solution = [];
    for k = 1:nCommodities
        Solution = [Solution solution{k}];
    end
    t = max(sum(Solution,2)./cap)

    alpha = 1/(iter+10);
    mu = mu + alpha * sum(Solution,2);

    H = 0.5 * eye(nArcs); f = - mu;
    mu = quadprog(H,f,[],[],cap',[1],zeros(nArcs,1),[]);
    iter = iter + 1;
end

end

