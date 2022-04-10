% relativeDirection_ctrlDistribution
%%
cd positionalData
load positionalDataByTrialType.mat
load platformLocations.mat
load frame.mat frameSize

%%
cd ..
cd physiologyData
cd direction

trialTypes = fieldnames(pos);

%  binSize = 5;
binSize = 15;
xAxis = 1:binSize:frameSize(1);
yAxis = 1:binSize:frameSize(2);
angleEdges = 0:15:360;

for tt = 1:length(trialTypes)       
    position = vertcat(pos.(trialTypes{tt})(:).dlc_XYsmooth);
    
    hd = [];
    plat = [];
    for t = 1:length(platformLocations.(trialTypes{tt}))
        
        hdTemp = vertcat(pos.(trialTypes{tt})(t).dlc_angle);
        hdNaN = find(isnan(hdTemp));
        hdIsN = find(~isnan(hdTemp));
        interpPhase = interpPhaseNEW(hdIsN, hdTemp(hdIsN), hdNaN);
        hdTemp(hdNaN) = interpPhase;
        hd = [hd; hdTemp];
        
        plat = [plat; vertcat(platformLocations.(trialTypes{tt}){t}.body)];
    end
    
    %%
    for p = 1:61 
        platInd = find(plat == p);
        if isempty(platInd)
            continue
        end
        posPlat = position(platInd,:);
        hdPlat = hd(platInd);
        
        posPlat(isnan(hdPlat),:) = [];
        hdPlat(isnan(hdPlat)) = [];
        
        purDirDistsTemp = histcounts(hdPlat, angleEdges);
        purDirDistsTemp = purDirDistsTemp./sum(purDirDistsTemp);
        purDirDists.(trialTypes{tt}){p} = purDirDistsTemp;    
        
        relDirDistsTemp = getRelDirDist(posPlat, hdPlat, xAxis, yAxis, angleEdges);
        relDirDists.(trialTypes{tt}){p} = relDirDistsTemp;
    end
end

save mrlFocus_ctrlDistribution_coarse.mat purDirDists relDirDists angleEdges
%%
function relDirDist = getRelDirDist(pos, hd, xAxis, yAxis, angleEdges)

xDistance = (xAxis - pos(:,1))';
xDistance = reshape(xDistance, [length(xAxis), 1, length(hd)]);
xDistance = repmat(xDistance, 1, length(yAxis), 1);

yDistance = (yAxis - pos(:,2))';
yDistance = reshape(yDistance, [1, length(yAxis), length(hd)]);
yDistance = repmat(yDistance, length(xAxis), 1, 1);

% euclidDist = sqrt(xDistance.^2 + yDistance.^2);
dir2goal = rad2deg(atan2(xDistance, yDistance));
dir2goal(dir2goal < 0) = dir2goal(dir2goal < 0) + 360;

spikeHD_XP = reshape(hd, 1, 1, length(hd));
spikeHD_XP = repmat(spikeHD_XP, length(xAxis), length(yAxis));

dirRel2Goal = spikeHD_XP - dir2goal;
dirRel2Goal(dirRel2Goal < 0) = 360 + dirRel2Goal(dirRel2Goal < 0);

relDirDist = NaN(length(xAxis), length(yAxis), length(angleEdges)-1);
for x = 1:length(xAxis)
    for y = 1:length(yAxis)
        distTemp = histcounts(dirRel2Goal(x,y,:), angleEdges);
        distTemp = distTemp./sum(distTemp);
        relDirDist(x,y,:) = distTemp;
    end
end
end

