function [fval,solution,dualCommodities, dualStations, dualArcs]=restricted_Master(pathsK_nodes, pathsK, nStations, node, nArcs, arc, nCommodities, commodities, cost)


number=0; % number of variables
% objective
for i=1:nCommodities
    number = number + size(pathsK{i},1);
end
c = zeros(number+1,1);
c(end) = 1;

% Flow conservation constraints:
count = 0;
Aeq = zeros(nCommodities,number+1); % linear constraint matrix
beq = zeros(nCommodities,1);
for i=1:nCommodities
    for j=1:size(pathsK{i},1)
        count = count + 1;
        Aeq(i,count) = 1;
    end
    beq(i) = commodities(i,3);
end

% capacity conservation constraints:
count = 0;
Aineq = [];
bineq = zeros(nArcs,1);

for j=1:nCommodities
%     for i=1:nArcs
%         for k=1:size(pathsK{j},1)
%             if pathsK{j}(k,i) > 0
%                 Aineq(i,count+k) = 1;
%             end
%         end
%     end
%     count = count + size(pathsK{j},1);
Aineq = [Aineq, pathsK{j}'];
end

Aineq = Aineq > 0;


Aineq = [Aineq -arc(:,4)];
% Aeq
% Aineq

lb = zeros(number+1,1);
ub = inf(number+1,1);
ub(end) = 1;

% solve
[solution,fval,~,~,lambda]=linprog(c,Aineq,bineq,Aeq,beq,lb,ub);
dualCommodities = lambda.eqlin;
dualStations = zeros(nStations,1);
dualArcs = lambda.ineqlin;

end



    