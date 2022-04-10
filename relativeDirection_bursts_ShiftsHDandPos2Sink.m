% relativeDirection_bursts_ShiftsHDandPos2Sink
%%
cd positionalData
load mazeCorners.mat
load goalCoordinates.mat
load platformLocations.mat
load positionalDataByTrialType.mat

cd ..

cd physiologyData

load burstSpikes.mat burstSpikes
% load burstSpikes_May2021.mat burstSpikes
load units.mat sample_rate


cd direction
load mrlFocus_bursts_coarse.mat
load mrlFocus_ctrlDistribution_coarse.mat angleEdges

trialTypes = fieldnames(pos);
%% get platform location for each frame - will use later to assign platform
% locations to individual spikes
% set min shift
minShift = 30 * 60; % in frames, 30 fps, so ~ 1 min

binSize = 15;
xAxis = 1:binSize:frameSize(1);
yAxis = 1:binSize:frameSize(2);
histBinCentres = deg2rad((angleEdges(2)-angleEdges(1))/2 + angleEdges(1:end-1));
nShuffles = [50, 100, 200, 500, 1000];

for tt = 1:length(trialTypes)
    samples = vertcat(pos.(trialTypes{tt})(:).sample);
    samples = round(samples);
    nFrames = length(samples);
    
    hd = [];
    plats = [];
    position = [];
    for t = 1:length(platformLocations.(trialTypes{tt}))
        plats = [plats; ...
            vertcat(platformLocations.(trialTypes{tt}){t}.body)];
        position = [position; ...
            pos.(trialTypes{tt})(t).dlc_XYsmooth];
        
        catAngleTemp = pos.(trialTypes{tt})(t).dlc_angle;
        nanAngles = find(isnan(catAngleTemp));
        ts = pos.(trialTypes{tt}).ts;
        if ~isempty(nanAngles)
            notNaN = find(~isnan(catAngleTemp));
            catAngleTemp(nanAngles) = rad2deg(interpPhaseNEW(ts(notNaN), ...
                deg2rad(catAngleTemp(notNaN)), ts(nanAngles)));
        end
        hd = [hd; catAngleTemp];
    end
    
    maxShift = nFrames - minShift;  
    
    shifts = [0; (randsample(maxShift - minShift, 1000) + minShift)];    
    %% get list of units to use
    nCells = length(signifCellsNorm.(trialTypes{tt}));
    mrlTemp = NaN(1001, nCells);
    directionsTemp = NaN(1001, nCells);
    coors = NaN(nCells, 2);
    
    spikeSamples = cell(nCells, 1);
    spikeFrames = cell(nCells, 1); 
    % spikeHD = cell(nCells, 1);
    nSpikes = cell(nCells,1);
    for c = 1:nCells
        unitName = signifCellsNorm.(trialTypes{tt}){c};
        
        mrlInd = find(strcmp(unitName, {mrlFocus(:).unit}));        
        coors(c,:) = mrlFocus(mrlInd).(trialTypes{tt}).norm.coor;        
        
        spikeSamples{c} = vertcat(burstSpikes.(unitName).(trialTypes{tt})(:).samples);
        
        % spikeHD{c} = NaN(length(spikeSamples{c}), 1);
        nSpikes{c} = NaN(length(spikeSamples{c}), 1);
        spikeFrames{c} = cell(length(spikeSamples{c}), 1);
        for b = 1:length(spikeSamples{c})
            spikeFrames{c}{b} = round(interp1(samples, ...
                1:length(samples), double(spikeSamples{c}{b})));
            % burstHD = hd(spikeFrames{c}{b});
            % hdTemp = rad2deg(circ_mean(deg2rad(burstHD)));
            % if hdTemp < 0
            %     hdTemp = hdTemp +360;
            % end
            % spikeHD{c}(b) = hdTemp;
            nSpikes{c}(b) = length(spikeSamples{c}{b});
        end
    end    
    
    signifCells_Shift.(trialTypes{tt}) = {};
    %%
    parfor s = 1:length(shifts)
        if s == 1
            platsShifted = plats;
            posShifted = position;
            hdShifted = hd;
        else
            platsShifted = plats([shifts(s):end, 1:shifts(s)-1]);
            posShifted = position([shifts(s):end, 1:shifts(s)-1], :);
            hdShifted = hd([shifts(s):end, 1:shifts(s)-1]);
        end
        
        relDirDists = cell(61,1);
        for p = 1:61
            platInd = find(platsShifted == p);
            if isempty(platInd)
                continue
            end
            posPlat = posShifted(platInd,:);
            hdPlat = hdShifted(platInd);
            
            posPlat(isnan(hdPlat),:) = [];
            hdPlat(isnan(hdPlat)) = [];
            
            relDirDists{p} = getRelDirDist(posPlat, hdPlat, xAxis, yAxis, angleEdges);
        end
        %%        
        for c = 1:nCells
            
            spikePos = [];
            spikeHD = [];
            for b = 1:length(spikeFrames{c})
                
                burstPos = posShifted(spikeFrames{c}{b}, :);
                burstPos = mean(burstPos, 1);
                
                spikePos = [spikePos; repmat(burstPos, nSpikes{c}(b), 1)];
                
                hdTemp = hdShifted(spikeFrames{c}{b});
                hdTemp = rad2deg(circ_mean(deg2rad(hdTemp)));
                if hdTemp < 0
                    hdTemp = hdTemp +360;
                end                
                spikeHD = [spikeHD; repmat(hdTemp, nSpikes{c}(b), 1)];
            end
            
            spikeFramesTemp = vertcat(spikeFrames{c}{:});
                     
            %%
            xInd = find(xAxis == coors(c,1));
            yInd = find(yAxis == coors(c,2));
            
            relDirDistsTemp = cell(61, 1);
            for p = 1:length(relDirDistsTemp)
                if isempty(relDirDists{p})
                    continue
                end
                relDirDistsTemp{p} = squeeze(relDirDists{p}(xInd, yInd, :))';
            end
            
            totalDist = ctrlDistributions(spikeFramesTemp, platsShifted, ...
                relDirDistsTemp);
            
            counts = dirHistCounts(spikePos, ...
                spikeHD, coors(c,1), coors(c,2), angleEdges, histBinCentres);
            
            normDist = counts./totalDist;
            sumNormDist = sum(normDist);
            normDistFactor = length(spikeHD) ./ sumNormDist;
            normDist = normDist .* normDistFactor;
            mrlDataNorm = mrlDir(normDist, histBinCentres);
            
            mrlTemp(s, c) = mrlDataNorm.mrl;
            directionsTemp(s, c) = mrlDataNorm.dir;
        end      
    end
    mrls.(trialTypes{tt}) = mrlTemp;
    directions.(trialTypes{tt}) = directionsTemp;
    shiftsUsed.(trialTypes{tt}) = shifts;
    %% put data into mrl structure
    for c = 1:nCells
        mrlReal = mrls.(trialTypes{tt})(1, c);
        
        mrlShift = mrls.(trialTypes{tt})(2:end, c);
        mrlShift(isnan(mrlShift)) = [];
        mrlShift = sort(mrlShift, 'descend');
        nShifts = length(mrlShift);
        sigLevel5 = floor(nShifts * 0.05);
        
        signifCells_Shift.(trialTypes{tt})(c).unit = ...
            signifCellsNorm.(trialTypes{tt}){c};
        
        signifCells_Shift.(trialTypes{tt})(c).mrl = mrlReal;
        
        if mrlReal > mrlShift(sigLevel5)
            signifCells_Shift.(trialTypes{tt})(c).significant = 'yes';
            
        else
            signifCells_Shift.(trialTypes{tt})(c).significant = 'no';
        end
        
        signifCells_Shift.(trialTypes{tt})(c).CI95 = ...
            mrlShift(sigLevel5);
        
        sigLevel2P5 = floor(nShifts * 0.025);
        signifCells_Shift.(trialTypes{tt})(c).CI97P5 = ...
            mrlShift(sigLevel2P5);
        
        if length(mrlShift) == 1000
            sigLevelZZ1 = floor(nShifts * 0.001);
            signifCells_Shift.(trialTypes{tt})(c).CI99P9 = ...
                mrlShift(sigLevelZZ1);
        end
        
        signifCells_Shift.(trialTypes{tt})(c).dir = directions.(trialTypes{tt})(1, c);
    end
end
%%
save mrlFocus_bursts_ShiftsHDandPos2Sink.mat angleEdges histBinCentres ...
    signifCells_Shift mrls coors directions shiftsUsed

%%
function relDirDist = getRelDirDist(pos, hd, xAxis, yAxis, angleEdges)

xDistance = (xAxis - pos(:,1))';
xDistance = reshape(xDistance, [length(xAxis), 1, length(hd)]);
xDistance = repmat(xDistance, 1, length(yAxis), 1);

yDistance = (yAxis - pos(:,2))';
yDistance = reshape(yDistance, [1, length(yAxis), length(hd)]);
yDistance = repmat(yDistance, length(xAxis), 1, 1);

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
%%
function totalDist = ctrlDistributions(spikeFrameInd, plat, dirDists)

spikePlats = plat(spikeFrameInd);
totalDist = [];

for p = 1:length(dirDists)
    nSpikesPerPlatform = length(find(spikePlats == p));
    if nSpikesPerPlatform == 0
        continue
    end
    
    dist = dirDists{p};
    dist = dist * nSpikesPerPlatform;
    
    if isempty(totalDist)
        totalDist = dist;
    else
        totalDist = totalDist + dist;
    end
end
end
%%
function histCounts = dirHistCounts(pos, hd, xAxis, yAxis, angleEdges, ...
    histBinCentres)

xDistance = (xAxis - pos(:,1))';
yDistance = (yAxis - pos(:,2))';

if length(xAxis) > 1
    xDistance = reshape(xDistance, [length(xAxis), 1, length(hd)]);
    xDistance = repmat(xDistance, 1, length(yAxis), 1);
    
    yDistance = reshape(yDistance, [1, length(yAxis), length(hd)]);
    yDistance = repmat(yDistance, length(xAxis), 1, 1);
end

dir2goal = rad2deg(atan2(xDistance, yDistance));
dir2goal(dir2goal < 0) = dir2goal(dir2goal < 0) + 360;

if length(xAxis) == 1
    dirRel2Goal = hd - dir2goal';
    dirRel2Goal(dirRel2Goal < 0) = 360 + dirRel2Goal(dirRel2Goal < 0);
    histCounts = histcounts(dirRel2Goal, angleEdges);
else
    spikeHD_XP = reshape(hd, 1, 1, length(hd));
    spikeHD_XP = repmat(spikeHD_XP, length(xAxis), length(yAxis));
    dirRel2Goal = spikeHD_XP - dir2goal;
    dirRel2Goal = permute(dirRel2Goal, [3, 2, 1]);
    
    dirRel2Goal(dirRel2Goal < 0) = 360 + dirRel2Goal(dirRel2Goal < 0);
    
    histCounts = NaN(length(histBinCentres), ...
        size(dirRel2Goal, 2), size(dirRel2Goal, 3));
    for x = 1:length(xAxis)
        for y = 1:length(yAxis)
            histCounts(:, y, x) = ...
                histcounts(dirRel2Goal(:, y, x), angleEdges);
        end
    end
end
end
%%
function mrlData = mrlDir(histCounts, histBinCentres)

mrlData.mrl = circ_r(histBinCentres', histCounts');
mrlData.dir = circ_mean(histBinCentres', histCounts');
end






