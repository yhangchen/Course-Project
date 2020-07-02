problem='assad1.5k';

[nStations, nodes, nArcs, arc, nCommodities, commodities, cost]=readdata(problem);

lambda = ones(nArcs,1);

% lagrange relaxiation
pathsK_nodes = cell(nCommodities,1);
pathsK = cell(nCommodities,1);
solution = cell(nCommodities,1);
for iter=1:1000
    iter
    alpha=1/i;
    beta=1/i;
    % for each k, compute the lagrange primal problem
    total_cap = 0;
    for k=1:nCommodities
        cost_k = -1*ones(nStations, nStations, 2);
        for i=1:nArcs
            cost_k(arc(i,1), arc(i,2),:) = [arc(i,3)*lambda(i), i];
        end
        capacity_node=[Inf];
%        [pathsK_nodes{k}, pathsK{k}, solution{k}]=columnGeneration(nStations, nodes, nArcs, arc, nCommodities, commodities, cost_k, capacity_node);
        [pathsK_nodes(k), pathsK(k), solution{k}]=columnGeneration(nStations, nodes, nArcs, arc, 1, commodities(k,:), cost_k, capacity_node);
    end
    % solve for t
    % calculate Rx
    Rx=zeros(nArcs,1);
    for l=1:nArcs
        sum=0; % R[l]x
        for k=1:nCommodities  %R_k[l]x_k
            pathnumber=size(pathsK{k},1);
            for j=1:pathnumber % count for all paths of commodity k.
                if pathsK{k}(j,l)~=0
                    sum=sum+solution{k}(j)*commodities(k,3);
                end
            end
        end
        Rx(l)=sum;
    end
    
    %   calculate c_l
    cap_arc=arc(:,4);
   
 
    % update lambda
    lambda=lambda+alpha*Rx;
    H = 0.5 * eye(size(lambda,1)); f = - lambda;
    mu = quadprog(H,f,[],[],cap_arc',[1],zeros(nArcs,1),[]);
    iter = iter + 1;
        t = max(Rx./cap_arc)

end
t
    
        
        
