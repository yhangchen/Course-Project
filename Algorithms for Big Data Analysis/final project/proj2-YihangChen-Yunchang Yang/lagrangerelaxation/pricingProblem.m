function [reducedCost, bestPath_nodes, forCommodity, bestPath_arcs]=pricingProblem(dualCommodities, dualStations, dualArcs,nStations, node, nArcs, arc, nCommodities, commodities, cost, capacity_node)

%     """ Return the paths corresponding to the most negative reduced cost of each commodity
% 
%         1. For each commodity, solve the pricing problem using Dijkstra's algorithm
%             1.1. If a path with a negative reduced cost is found, add it to the solution set
%         2. Return the solution set, and the corresponding commodities
%     """ 

% pathsK_nodes is a 1d cell array, the first dimension is the id of
% commodity,pathsK_nodes{i} is a matrix of size p*(nstation+1), p is the
% number of paths of commodity i, each row contains information about the
% nodes on the path and cost.

% bestPath_nodes is a 1d array, the length is the same as forCommodity. So
% is bestPath_arcs.

reducedCost = 0; % Contains a negative value if we found a path with negative reduced cost
forCommodity = []; % id of the commodity which will receive the new paths, if any
bestPath_nodes = {}; % Contains the paths which will be added (node_id's and cost)
bestPath_arcs = {}; % Contains the paths which will be added (arc_id's)

shortest_from_node = ones(nStations, nStations)*Inf; % Total cost of the shortest path from index
previous_nodes = ones(nStations, nStations)*-1; % Nodes id to go backward
previous_arcs = ones(nStations, nStations)*-1; % Arcs id to go backward


% -- Compute the new cost vector -- 
cost_for_commo = cost;
for i=1:nStations
    for j=1:nStations
        if cost_for_commo(i,j,1)~=0
            cost_for_commo(i,j,1) = cost_for_commo(i,j,1)+dualStations(j);
        end
    end
end

for i=1:length(dualArcs)
    if dualArcs~=0
        if cost_for_commo(arc(i,1),arc(i,2),1)>0
            cost_for_commo(arc(i,1),arc(i,2),1)=cost_for_commo(arc(i,1),arc(i,2),1)+dualArcs(i);
        end
    end
end

% -- Run Dijkstra's -- 
count = 1;
for i=1:nCommodities
    origin=commodities(i,1);
    dest=commodities(i,2);
    
    if ~isequal(previous_nodes(origin,:),ones(1,nStations)*-1) % Haven't run Dijkstra yet for this origin
        [shortest_from_node(origin,:), previous_nodes(origin,:), previous_arcs(origin,:)] = dijkstra(origin, cost_for_commo, [], []);
    end
    
    currentReducedCost = (shortest_from_node(origin, dest) * commodities(i,3)) + dualCommodities(i); %Compute the reducedCost
    
    if currentReducedCost<0
        % Construct the new path
%         newPath_nodes = zeros(1,nStations+1);
        newPath_arcs = zeros(1,nArcs);
        
        % By going backward from destination
        while dest~=origin
%             newPath_nodes(dest)=1 * commodities(i,3);
%             newPath_nodes(nStations+1)=newPath_nodes(nStations+1)+arc(previous_arcs(origin,dest),3); % Compute the cost on the fly
            
            newPath_arcs(previous_arcs(origin,dest)) = 1 * commodities(i,3);
            
            dest = previous_nodes(origin, dest);
        end
        
        forCommodity=[forCommodity i];
%         bestPath_nodes{count}=newPath_nodes;
        bestPath_arcs{count}=newPath_arcs;
        count = count + 1;
    end
end
        
        
    
    
