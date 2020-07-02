function [labels, precedent_nodes, precedent_arcs] = dijkstra(origin, cost, mostUsedNode, mostUsedArc, nStation, arc)
% % Disconnect the mostUsedNode from all others
% if ~isempty(mostUsedNode)
%     for jj = 1:length(mostUsedNode)
%         for ii = 1:nStation
%             cost(mostUsedNode(jj),ii,1) = -1;
%             cost(mostUsedNode(jj),ii,2) = -1;
%         end
%     end
% end
% Disconnect the mostUsedArc from all others 
if ~isempty(mostUsedArc)
    for jj = 1:length(mostUsedArc)
        arcjj = arc(mostUsedArc(jj),:);
        cost(arcjj(1),arcjj(2),1) = -1;
        cost(arcjj(1),arcjj(2),2) = -1;
    end
end

labels = inf(nStation,1);
labels(origin) = 0;

precedent_nodes = -ones(1,nStation);
precedent_arcs = -ones(1,nStation);

toTreat = [origin];

while true
    [~,currentIndex] = min(labels(toTreat));
    currentNode = toTreat(currentIndex);
    neighbors = find(cost(currentNode,: ,2)>-1);
    for ii = 1:length(neighbors)
        neighbor = neighbors(ii);
        if labels(neighbor) > labels(currentNode) + cost(currentNode, neighbor,1)
            labels(neighbor) = labels(currentNode) + cost(currentNode, neighbor,1);
            precedent_nodes(neighbor) = currentNode;
            precedent_arcs(neighbor) = cost(currentNode, neighbor, 2);
            toTreat(end+1) = neighbor;
        end
    end
    toTreat(currentIndex) = [];
    if isempty(toTreat)
        break
    end
end
end



