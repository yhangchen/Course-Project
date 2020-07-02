function baseline(problem)
[nStations, nodes, nArcs, arc, nCommodities, commodities, cost]=readdata(problem);

B = zeros(nStations,nArcs);
for i = 1:nArcs
    B(arc(i,1),i) = 1;
    B(arc(i,2),i) = -1;
end

b = cell(1,nCommodities);
b_full = []; B_full = []; subB_full = [];
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
    b_full = [b_full; btmp];
    B_full = blkdiag(B_full,B);
end
cap = arc(:,4);

B_full(:,end+1) = 0;
Aineq = repmat(eye(nArcs),1,nCommodities);
Aineq(:,end+1) = -cap;
bineq = zeros(size(Aineq,1),1);
lb = zeros(size(B_full,2),1);
c = zeros(size(B_full,2),1);
c(end) = 1;

ub = inf(size(B_full,2),1);
ub(end) = 1;

[solution,objfunction,~,~,~]=linprog(c,Aineq,bineq,B_full,b_full,lb,ub);

objfunction