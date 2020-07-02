function [pathsK_nodes, pathsK, solution]=columnGeneration(nStations, node, nArcs, arc, nCommodities, commodities, cost, capacity_node)
% pathsK_nodes[k,i,j] with k = commodity_id // i = path_id of commodity k // j = node_id of the path 
% pathsK_nodes is a 1d cell array, the first dimension is the id of
% commodity,pathsK_nodes{i} is a matrix of size p*(nstation+1), p is the
% number of paths of commodity i, each row contains information about the
% nodes on the path and cost.
pathsK_nodes = cell(nCommodities,1);

% pathsK[k,i,j] with k = commodity_id // i = path_id of commodity k // j = arc_id of the path 
pathsK = cell(nCommodities,1); % Contains all initial paths for each commodity (by arcs)

[pathsK_nodes, pathsK] = getInitSet(nStations, node, nArcs, arc, nCommodities, commodities, cost, capacity_node); % Get the initial set of variable through Dijkstra and "route through best possible shortest path" method

%t_time_init_set = time.time() # Time - 2 Start optimizing the problem

while 1
    [obj_Function, solution, dualCommodities, dualStations, dualArcs] = restricted_Master(pathsK_nodes, pathsK, nStations, node, nArcs, arc, nCommodities, commodities, cost, capacity_node); % Solve the restricted master problem
    [reducedCost, newPath, forCommodity, newPath_arc] = pricingProblem(dualCommodities, dualStations, dualArcs,nStations, node, nArcs, arc, nCommodities, commodities, cost, capacity_node); % Solve the pricing problem
    if(abs(reducedCost) < 10^-21)
        
        % return solution
        return 

    else
        for j=1:length(newPath_arc)
%             pathsK_nodes(forCommodity(j))=[pathsK_nodes(forCommodity(j));newPath(j)];
            pathsK(forCommodity(j))=[pathsK(forCommodity(j));newPath_arc(j)];
        end
    end
end

        