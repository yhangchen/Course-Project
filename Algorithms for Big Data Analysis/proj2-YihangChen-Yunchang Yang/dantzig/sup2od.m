function suppath = sup2od(filename, nStations)
suppath = load([filename,'.sup']);
addpath = []; not_od_ind = [];
for i = 1:size(suppath,1)
    if suppath(i,1) < 0 || suppath(i,2) < 0
        not_od = suppath(i,:);
        not_od_ind = [not_od_ind i];
        if not_od(1) < 0 && not_od(2) > 0
            ii = i - 1;
            preset = [];
            while ii > 0 && suppath(ii,2) == suppath(i,2)
                preset = [preset suppath(ii,1)];
                ii = ii - 1;
            end
              
            for j = 1:nStations
                if ~ismember(j, preset)
                    addpath = [addpath; j not_od(2) not_od(3) not_od(4)];
                end
            end
        elseif not_od(1) > 0 && not_od(2) < 0
            ii = i - 1;
            preset = [];
            while ii > 0 && suppath(ii,1) == suppath(i,1)
                
                preset = [preset suppath(ii,2)];
                ii = ii - 1;
            end
              
            for j = 1:nStations
                if ~ismember(j, preset)
                    addpath = [addpath; not_od(1) j not_od(3) not_od(4)];
                end
            end
        end
    end
end

suppath(not_od_ind,:) = [];
% suppath = [suppath;addpath];

end