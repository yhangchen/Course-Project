function [fval,solution,dualCommodities, dualStations, dualArcs]=restricted_Master(pathsK_nodes, pathsK, nStations, node, nArcs, arc, nCommodities, commodities, cost, capacity_node)

number=0; % number of variables
% objective
c=[];
for i=1:nCommodities
    for j=1:size(pathsK{i},1)
        cc = sum(pathsK{i}(j,:) ~= 0);
        c=[c cc*commodities(i,3)];
%         c=[c pathsK_nodes{i}(j,nStations+1)];
        number = number +1;
    end
end
        

% Flow conservation constraints :
count = 1;
Aeq = zeros(nCommodities,number); % linear constraint matrix
beq = zeros(nCommodities,1);
for i=1:nCommodities
    ind=[];
    val=[];
    for j=1:size(pathsK{i},1)
        Aeq(i,count) = 1;
        count = count + 1;
    end
    beq(i) = 1;
    
end



lb = zeros(number,1);

% solve
options = optimoptions('linprog','Display','none');
[solution,fval,~,~,lambda]=linprog(c,[],[],Aeq,beq,lb,[],options);
dualCommodities = lambda.eqlin;
dualStations = zeros(nStations,1);
dualArcs = zeros(nArcs,1);



    