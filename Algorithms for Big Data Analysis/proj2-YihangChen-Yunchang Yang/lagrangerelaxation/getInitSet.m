function [pathsK_nodes,pathsK]=getInitSet(nStations, node, nArcs, arc, nCommodities, commodities, cost, capacity_node)

%    #   -------- GET INITIAL SET OF PROMISSING VARIABLES --------

shortest_paths_original = zeros(nStations, nStations); % Represent the shortest path distance between each nodes in the original network
previous_nodes_original = zeros(nStations, nStations); % Represent the previous node of each node of the shortest path from index
previous_arcs_original = zeros(nStations, nStations); % Represent the previous arc of each node of the shortest path from index

distance_done = []; % Nodes for which dijkstra's algorithm has already been run
sorting = 0; % Number of iteration over all commodities
    
outdegree = zeros(nCommodities,1); % Contain the outdegree of the shortest path of each commodity
sortedCommo = 1:nCommodities; % Contains the indices by which we are going to iterate over the commodities

maxIter = 1000;

while sorting < maxIter
    current_capacities_arcs = arc(:,4); % To verify the capacities of the arcs are not violated
    pathsK_nodes=cell(1,nCommodities);
    pathsK=cell(1,nCommodities);
    
    unfeasible = 0;
    count=0;
    
    not_add = 1;
    
    for ii=1:length(sortedCommo)
        i = sortedCommo(ii);
        count = count + 1;
        current_cost=cost;
        to_Ignore_node = [];
        to_Ignore_arc = [];
        stopping = 0;
        commodities_in = commodities;
        while ~stopping
            start=commodities(i,1);
            dest=commodities(i,2);
            if sorting == 0 % 1st iteration -> Compute the shortest path
                if ~ismember(start,distance_done)
                    [shortest_path_original(start,:),previous_nodes_original(start,:), previous_arcs_original(start,:)] = dijkstra(start, current_cost, to_Ignore_node, to_Ignore_arc, nStations, arc);
                    distance_done=[distance_done start];
                end
                
                % Set the current shortest path values  
                shortest_paths= shortest_paths_original(start,:);
                previous_nodes=previous_nodes_original(start,:);
                previous_arcs= previous_arcs_original(start,:);
            elseif not_add==1 % # A node/arc has to be ignored, re-run dijkstra's algorithm
                [shortest_paths, previous_nodes, previous_arcs] = dijkstra(start, current_cost, to_Ignore_node, to_Ignore_arc, nStations, arc);
            else %# If no nodes/arcs have to be ignored, take the original shortest path value
                shortest_paths= shortest_paths_original(start,:);
                previous_nodes=previous_nodes_original(start,:);
                previous_arcs= previous_arcs_original(start,:);
            end
            
            
            if shortest_paths(dest) < Inf % Check the feasibilty of the commodity
                to_Ignore_arc_in=[];
                
                % -- Generate the new path --
                newPath_arcs=zeros(1,nArcs);
                allocation_ratio = 1;
                encode_allocation_ratio = [];
                while dest ~= start
                    if sorting>0
                        if current_capacities_arcs(previous_arcs(dest)) - commodities_in(i,3) < 0
                            %fprintf(['block at arc : \n', num2str(previous_arcs(dest))]);
                            allocation_ratio = min(allocation_ratio,current_capacities_arcs(previous_arcs(dest))/commodities_in(i,3));
                            encode_allocation_ratio = [encode_allocation_ratio, allocation_ratio];
%                             if ~ismember(previous_arcs(dest), to_Ignore_arc_in)
                            to_Ignore_arc_in=[to_Ignore_arc_in;previous_arcs(dest)];
%                             end
                        end
                    end

                    
                    if sorting==0
                        if previous_nodes(dest)>0
                            outdegree(i) = outdegree(i) + node(previous_nodes(dest),2);
                        end
                    end
                    newPath_arcs(previous_arcs(dest)) = commodities_in(i,3);
                    
                    dest = previous_nodes(dest);
                   
                end
                if sorting > 0
                    to_Ignore_ind = (encode_allocation_ratio-allocation_ratio) < eps;
                    to_Ignore_arc_in = to_Ignore_arc_in(to_Ignore_ind);
                    commodities_in(i,3) =  (1-allocation_ratio)*commodities_in(i,3);
                    newPath_arcs = allocation_ratio * newPath_arcs;
                
                    if allocation_ratio
                        pathsK{i}=[pathsK{i};newPath_arcs];
                    end
                    current_capacities_arcs = current_capacities_arcs - newPath_arcs';
                end
            if allocation_ratio == 1
                break
            else
                to_Ignore_arc = [to_Ignore_arc;to_Ignore_arc_in];
            end

            else    %# There is no feasible path between the origin and destination
                unfeasible = 1;
                argwhere = sortedCommo==i;
                sortedCommo(argwhere)=[];
                sortedCommo = [i;sortedCommo];
%                 fprintf('unfeasible, try again\n');
                break
            end
        end
        
    end

        
    if sorting==0
        [~,sortedCommo] = sort(outdegree);
    end
%     sortedCommo = sortedCommo(randperm(length(sortedCommo)));
    sorting = sorting+1;
    if unfeasible == 0 && sorting>1
        return 
    end
end

if sorting == maxIter
    error('infeasible')
end

end
                         
                    
                            
                
    

