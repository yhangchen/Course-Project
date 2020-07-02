problem='assad3.7k';
t0 = clock;
[nStations, nodes, nArcs, arc, nCommodities, commodities, cost]=readdata(problem);
t1 = clock;
solution = cell(nCommodities,1);
[pathsK_nodes, pathsK] = getInitSet(nStations, nodes, nArcs, arc, nCommodities, commodities, cost, []);
t2 = clock;

while true
    [obj_Function, solution, dualCommodities, dualStations, dualArcs] = restricted_Master(pathsK_nodes, pathsK, nStations, nodes, nArcs, arc, nCommodities, commodities, cost);
    
    [reducedCost, newPath, forCommodity, newPath_arc] = pricingProblem(dualCommodities, dualStations, dualArcs, nStations, nodes, nArcs, arc, nCommodities, commodities, cost);
    if abs(reducedCost) < eps
        t3 = clock;
       % check constraint 2
       dualArcs'*arc(:,4) >= 1-1e-6

        fprintf("\nt = " + string(obj_Function) + '\n')
        fprintf("\nLoading data duration : "+ etime(t1,t0)+ '\n')
        fprintf("Finding initial set duration : "+ etime(t2,t1)+ '\n')
        fprintf("Solving problem duration : "+ etime(t3,t2)+ '\n')
        fprintf("\nTotal duration : "+ etime(t3,t0)+ '\n')
        break
    
    else
        for j = 1:length(newPath)
            pathsK{forCommodity(j)} = [pathsK{forCommodity(j)};newPath_arc{j}];
        end
    end

end
obj_Function