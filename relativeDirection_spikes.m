% relativeDirection_spikes
%%
cd positionalData
load mazeCorners.mat
load goalCoordinates.mat
load platformLocations.mat
load positionalDataByTrialType.mat

cd ..

cd physiologyData

cd placeFieldData
load rateMaps.mat rateMaps frameSize
cd ..

cd direction
load mrlFocus_ctrlDistribution_coarse.mat

trialTypes = fieldnames(pos);
%% get platform location for each frame - will use later to assign platform 
% locations to individual spikes
samples = cell(1,2);
plat = cell(1,2);
for tt = 1:length(trialTypes)
    samples{tt} = vertcat(pos.(trialTypes{tt})(:).sample);
    samples{tt} = round(samples{tt});

    for t = 1:length(platformLocations.(trialTypes{tt}))       
        plat{tt} = [plat{tt}; vertcat(platformLocations.(trialTypes{tt}){t}.body)];
    end
end
%%
nCells = length(rateMaps);
cellCounter = 0;

histBinCentres = deg2rad((angleEdges(2)-angleEdges(1))/2 + angleEdges(1:end-1));

% binSize = 5;
binSize = 15;
xAxis = 1:binSize:frameSize(1);
yAxis = 1:binSize:frameSize(2);

nShuffles = [50, 100, 200, 500, 1000];
 
for tt = 1:length(trialTypes)
    signifCellsNorm.(trialTypes{tt}) = {};
end

for c = 1:nCells
    if ~strcmp(rateMaps(c).neuronType, 'pyramid') 
        continue
    end
      
    cellFlag = 0;
    
    unitName = rateMaps(c).unit;
    
    for tt = 1:length(trialTypes)
        spikePos = rateMaps(c).(trialTypes{tt}).spikePos;
        spikeHD = rateMaps(c).(trialTypes{tt}).spikeHD;
        spikeSamples = double(rateMaps(c).(trialTypes{tt}).spikeSamples);
        
        if length(spikeHD) < 500
            continue
        else
            cellFlag = cellFlag + 1;
        end
        
        if cellFlag == 1
            cellCounter = cellCounter + 1;
            mrlFocus(cellCounter).unit = unitName;
        end       
        %%
        platInd = interp1(samples{tt}, 1:length(samples{tt}), spikeSamples); clear spikeSamples
        platInd = round(platInd);
        spikePlats = plat{tt}(platInd); 
        
        platDist_Rel = [];
        totalDist_Rel = [];
        
        platDist_Pure = [];
        totalDist_Pure = [];
        
        for p = 1:length(relDirDists.(trialTypes{tt}))
            nSpikesPerPlatform = length(find(spikePlats == p));
            if nSpikesPerPlatform == 0
                continue
            end
            
            platDist_Rel = relDirDists.(trialTypes{tt}){p};            
            platDist_Rel = platDist_Rel * nSpikesPerPlatform;
            
            platDist_Pure = purDirDists.(trialTypes{tt}){p};            
            platDist_Pure = platDist_Pure * nSpikesPerPlatform;
            
            if isempty(totalDist_Rel)
                totalDist_Rel = platDist_Rel;
                totalDist_Pure = platDist_Pure;
            else
                totalDist_Rel = totalDist_Rel + platDist_Rel;
                totalDist_Pure = totalDist_Pure + platDist_Pure;
            end
        end
        
        totalDist_permute = permute(totalDist_Rel, [3, 2, 1]);
        %%
        dirRel2Goal_histCounts = dirHistCounts(spikePos, spikeHD, xAxis, yAxis, ...
            angleEdges, histBinCentres);
                
        normDist = dirRel2Goal_histCounts./totalDist_permute;
        sumNormDist = sum(normDist, 1);
        normDistFactor = length(spikeHD) ./ sumNormDist;
        normDistFactor = repmat(normDistFactor, length(histBinCentres), 1, 1);
        normDist = normDist .* normDistFactor;        
        mrlDataNorm = mrlRelDir(normDist, xAxis, yAxis, histBinCentres);  
        
        mrlFocus(cellCounter).(trialTypes{tt}).norm = mrlDataNorm;        
        %%
        hdCountsPure = histcounts(spikeHD, angleEdges);
        
        hdCountsNorm = hdCountsPure./totalDist_Pure;
        hdCountsNorm = hdCountsNorm./sum(hdCountsNorm);
        hdCountsNorm = hdCountsNorm * length(spikeHD);
        
        mrl = circ_r(histBinCentres', hdCountsNorm');
        direction = rad2deg(circ_mean(histBinCentres', hdCountsNorm'));
        if direction <  0
            direction = direction + 360;
        end
       [pval, z] = circ_rtest(histBinCentres', hdCountsNorm');
        
        mrlFocus(cellCounter).(trialTypes{tt}).pureNorm.mrl = mrl;
        mrlFocus(cellCounter).(trialTypes{tt}).pureNorm.dir = direction;
        mrlFocus(cellCounter).(trialTypes{tt}).pureNorm.distribution = ...
            hdCountsNorm;
        mrlFocus(cellCounter).(trialTypes{tt}).pureNorm.pval = pval;
        mrlFocus(cellCounter).(trialTypes{tt}).pureNorm.z = z;           
        %% shuffles
        mrlVal_norm = [];
        for sh = 1:length(nShuffles)
            nShufflesTemp = nShuffles(sh) - length(mrlVal_norm);
            mrlVal_Temp = NaN(nShufflesTemp, 1);
            
            if nShufflesTemp < 200
                sigLevel = floor(200/20);
            elseif nShufflesTemp == 200
                sigLevel = 20;
            else
                sigLevel = floor(1000/20);
            end
            
            parfor s = 1:nShufflesTemp
                spikeHDtemp = spikeHD;
                spikeHDtemp = randsample(spikeHDtemp, length(spikeHDtemp));
                
                dirRel2Goal_histCounts = [];
                
                dirRel2Goal_histCounts = dirHistCounts(spikePos, spikeHDtemp, ...
                    xAxis, yAxis, angleEdges, histBinCentres);
                
                normDist = dirRel2Goal_histCounts./totalDist_permute;
                sumNormDist = sum(normDist, 1);
                normDistFactor = length(spikeHD) ./ sumNormDist;
                normDistFactor = repmat(normDistFactor, length(histBinCentres), 1, 1);
                normDist = normDist .* normDistFactor;
                mrlNormShuffle = mrlRelDir(normDist, xAxis, yAxis, histBinCentres);
                mrlVal_Temp(s) = mrlNormShuffle.mrl;
            end
            
            mrlVal_norm = [mrlVal_norm; mrlVal_Temp];
            mrlVal_norm = sort(mrlVal_norm, 'descend');
            
            if mrlDataNorm.mrl < mrlVal_norm(sigLevel)
                sigLevel = floor(length(mrlVal_norm)/20);
                mrlFocus(cellCounter).(trialTypes{tt}).norm.CI95 = mrlVal_norm(sigLevel);
                break
            end
            
            if nShuffles(sh) == 1000
                mrlFocus(cellCounter).(trialTypes{tt}).norm.CI95 = mrlVal_norm(50);
                mrlFocus(cellCounter).(trialTypes{tt}).norm.CI97P5 = mrlVal_norm(25);
                mrlFocus(cellCounter).(trialTypes{tt}).norm.CI99P9 = mrlVal_norm(1);
                
                if mrlDataNorm.mrl > mrlFocus(cellCounter).(trialTypes{tt}).norm.CI95
                    signifCellsNorm.(trialTypes{tt}) = [signifCellsNorm.(trialTypes{tt}); ...
                        unitName];
                end
            end
        end
    end
end
%%
save mrlFocus_spikes_coarse.mat mrlFocus angleEdges histBinCentres signifCellsNorm     
% save mrlFocus_spikes_interneurons.mat mrlFocus angleEdges histBinCentres signifCellsNorm     
%%
for tt = 1:2
    signifCellsNorm.(trialTypes{tt}) = {};
    for c = 1:length(mrlFocus)
        if isempty(mrlFocus(c).(trialTypes{tt}))
            continue
        end
               
        unitName = {mrlFocus(c).unit};
        mrlDataNorm = mrlFocus(c).(trialTypes{tt}).norm;
        
        if mrlDataNorm.mrl > mrlDataNorm.CI95
            signifCellsNorm.(trialTypes{tt}) = [signifCellsNorm.(trialTypes{tt}); ...
                unitName];
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
function mrlData = mrlRelDir(histCounts, xAxis, yAxis, histBinCentres)

histBinCenRep = repmat(histBinCentres', 1, length(yAxis), length(xAxis));
mrl = circ_r(histBinCenRep, histCounts);
mrl = squeeze(mrl);
mrl_Lin = mrl(:);
[mrl_Max, mrlInd] = max(mrl_Lin);
[mrl_MaxCoor(1), mrl_MaxCoor(2)] = ...
    ind2sub(size(mrl), mrlInd);
mrlData.mrl = mrl_Max;

direction = circ_mean(histBinCenRep, histCounts);
direction = squeeze(direction);
mrlData.dir = rad2deg(direction(mrl_MaxCoor(1), mrl_MaxCoor(2)));
[pval, z] = circ_rtest(histBinCentres', histCounts(:, mrl_MaxCoor(1), mrl_MaxCoor(2)));
mrlData.distribution = squeeze(histCounts(:, mrl_MaxCoor(1), mrl_MaxCoor(2)))';
mrl_MaxCoor(1) = yAxis(mrl_MaxCoor(1));
mrl_MaxCoor(2) = xAxis(mrl_MaxCoor(2));
mrlData.coor = fliplr(mrl_MaxCoor);
mrlData.pval = pval;
mrlData.z = z;

end







