function [nStations, nodes, nArcs, arc, nCommodities, commodities, cost]=readdata(problem)

% nStations: number of nodes
% nodes: matrix nStations x 2, each row [index, outdegree];
% 
% nArcs: number of links (all bundled)
% arc: matrix nArcs x 4, each row [FromNode, ToNode, Cost(set to 1), Capacities]
% 
% nCommodities: number of OD pair
% commodities: matrix of nCommodities x 3, each row [FromNode, ToNode, bandwidth]
% 
% cost: matrix nStations x nStations x 2
% 
% cost(arc(i,1), arc(i,2),1) = arc(i,3); encode the cost
% cost(arc(i,1), arc(i,2),2) = i; index



nodepath=[problem,'.nod'];

node=load(nodepath);
nProducts = node(1);
nStations = node(2);
nArcs = node(3);
nBundles = node(4);

if nArcs ~= nBundles
    error('All links must be bundled link');
end

mutpath = [problem,'.mut'];
mut = load(mutpath);

arcpath = [problem,'.arc'];
initarc = load(arcpath); 
arc = [];
counter = zeros(nStations, nStations);
for i=1:size(initarc,1)
    row=initarc(i,:);
    if counter(row(1),row(2))==0
        if row(length(row))
            arc = [arc; row(1) row(2) 1 mut(row(end),2)];    %FromNode - ToNode - Cost - Capacities
        else
             error('All links must be bundled link');
        end
        counter(row(1),row(2))=1;
    end
end

if size(arc,1)~=nArcs
    fprintf('arc size error');
end

nodes = zeros(nStations,2);

for i=1:nStations
    nodes(i,1) = i;
end

% getcost
cost = -ones(nStations, nStations, 2);
% # Fill the cost variable by reading the arcs costs and indices

for i=1:nArcs
    cost(arc(i,1), arc(i,2),1) = arc(i,3);
    cost(arc(i,1), arc(i,2),2) = i;
end

%# def getOD(filepath):

try
suppath = [problem,'.od'];
commodities = load(suppath);
catch
 commodities = sup2od(problem, nStations);
end

nCommodities = size(commodities,1);

commodities(:,3)=[];

%    return nCommodities, commodities
%    #Compute the outdegree of the each node
%outdegree = np.count_nonzero(cost[:, :, 0], axis = 0)

outdegree=zeros(nStations,1);
for i=1:nStations
    outdegree(i)=nnz(counter(i,:));
    
%    Append it to the end of each entry
    
%     if(outdegree(i) == 0)
%         fprintf('NODE WITH INDEX : ', i ,'IS ISOLATED');
%     end
    nodes(i,2)=outdegree(i);

end
end



    



    