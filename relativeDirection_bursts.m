% relativeDirection_bursts
%%
cd positionalData
load mazeCorners.mat
load goalCoordinates.mat
load platformLocations.mat
load positionalDataByTrialType.mat

cd ..

cd physiologyData
load burstSpikes.mat burstSpikes

cd direction
%  load mrlFocus_ctrlDistribution.mat
load mrlFocus_ctrlDistribution_coarse.mat
%%
trialTypes = fieldnames(pos);

unitNames = fieldnames(burstSpikes);
nCells = length(unitNames);

cellCounter = 0;

histBinCentres = deg2rad((angleEdges(2)-angleEdges(1))/2 + angleEdges(1:end-1));

binSize = 15;
xAxis = 1:binSize:frameSize(1);
yAxis = 1:binSize:frameSize(2);

%nShuffles = 200;
nShuffles = [50, 100, 200, 500, 1000];

% signifCellsRaw.(trialTypes{1}) = {};
% signifCellsRaw.(trialTypes{2}) = {};

for tt = 1:length(trialTypes)
    signifCellsNorm.(trialTypes{tt}) = {};
end
for c = 1:nCells
    cellFlag = 0;
    
    for tt = 1:length(trialTypes)       
        bSpikesTemp = burstSpikes.(unitNames{c}).(trialTypes{tt});
        
        if isempty(vertcat(bSpikesTemp(:).samples)) || ...
                length(vertcat(bSpikesTemp(:).samples)) < 10
            continue
        end
        
        burstCounter = 0;
        nTrials = length(bSpikesTemp);
        burstAveragedSpike = [];
        
        for t = 1:nTrials
            if isempty(bSpikesTemp(t).samples)
                continue
            end
            
            nBurstsTemp = length(bSpikesTemp(t).samples);
            
            for b = 1:nBurstsTemp
                burstCounter = burstCounter + 1;
                
                if isnan(mean(bSpikesTemp(t).spikePos{b}, 1))
                    error
                end
                
                burstAveragedSpike(burstCounter,1).position = ...
                    mean(bSpikesTemp(t).spikePos{b}, 1);
                
                hdTemp = rad2deg(circ_mean(deg2rad(bSpikesTemp(t).spikeHD{b})));
                hdTemp(hdTemp < 0) = hdTemp(hdTemp < 0) + 360;
                burstAveragedSpike(burstCounter,1).hd = hdTemp;
                burstAveragedSpike(burstCounter,1).nSpikes = ...
                    length(bSpikesTemp(t).samples{b});
                
                samplesPlat = pos.(trialTypes{tt})(t).sample;
                samplesPlat = round(samplesPlat);
                plats = vertcat(platformLocations.(trialTypes{tt}){t}.body);
                burstSample = bSpikesTemp(t).samples{b}(round(length(bSpikesTemp(t).samples{b})/2));
                [~, platInd] = min(abs(double(burstSample)-samplesPlat));
                burstAveragedSpike(burstCounter,1).platform = plats(platInd);
            end
        end
        
        dataFields = {'position', 'hd', 'platform'};
        for f = 1:length(dataFields)
            burstUpSampled.(dataFields{f}) = [];
            for b = 1:length(burstAveragedSpike)                
                burstUpSampled.(dataFields{f}) = [burstUpSampled.(dataFields{f}); ...
                    repmat(burstAveragedSpike(b).(dataFields{f}), ...
                    burstAveragedSpike(b).nSpikes, 1)];
            end
        end
        
        if length(burstUpSampled.hd) < 500
            continue
        else
            cellFlag = cellFlag + 1;
        end
        
        if cellFlag == 1
            cellCounter = cellCounter + 1;
            mrlFocus(cellCounter).unit = unitNames{c};
            
            if cellCounter == 1 && tt == 2
                mrlFocus(cellCounter).hComb = [];
            end
        end
        %%
        platDist_Rel = [];
        totalDist_Rel = [];
        
        platDist_Pure = [];
        totalDist_Pure = [];
        
        for p = 1:length(relDirDists.(trialTypes{tt}))
            nSpikesPerPlatform = length(find(burstUpSampled.platform == p));
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
        dirRel2Goal_histCounts = dirHistCounts(burstUpSampled.position, ...
            burstUpSampled.hd, xAxis, yAxis, angleEdges, histBinCentres);
              
        normDist = dirRel2Goal_histCounts./totalDist_permute;
        sumNormDist = sum(normDist, 1);
        normDistFactor = length(burstUpSampled.hd) ./ sumNormDist;
        normDistFactor = repmat(normDistFactor, length(histBinCentres), 1, 1);
        normDist = normDist .* normDistFactor;
        mrlDataNorm = mrlRelDir(normDist, xAxis, yAxis, histBinCentres);
        mrlFocus(cellCounter).(trialTypes{tt}).norm = mrlDataNorm;
        %%
        hdCountsPure = histcounts(burstUpSampled.hd, angleEdges);
        
        hdCountsNorm = hdCountsPure./totalDist_Pure;
        hdCountsNorm = hdCountsNorm./sum(hdCountsNorm);
        hdCountsNorm = hdCountsNorm * length(burstUpSampled.hd);
        
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
                burstAvgTemp = burstAveragedSpike;
                spikeHDtemp = vertcat(burstAvgTemp(:).hd);
                spikeHDtemp = randsample(spikeHDtemp, length(spikeHDtemp));
                spikeHDtemp = num2cell(spikeHDtemp);
                [burstAvgTemp.hd] = spikeHDtemp{:};
                
                hdUp = [];
                posUp = [];
                
                for b = 1:length(burstAvgTemp)
                    hdUp = [hdUp; repmat(burstAvgTemp(b).hd, ...
                        burstAvgTemp(b).nSpikes, 1)];
                    posUp = [posUp; repmat(burstAvgTemp(b).position, ...
                        burstAvgTemp(b).nSpikes, 1)];
                end
                
                %             dirRel2Goal_histCounts_raw = dirHistCounts(posUp, ...
                %                 hdUp, coorX_raw, coorY_raw, angleEdges, histBinCentres);
                %             mrlVal_raw(s) = circ_r(histBinCentres', dirRel2Goal_histCounts_raw');
                %
                dirRel2Goal_histCounts_norm = dirHistCounts(posUp, ...
                    hdUp, xAxis, yAxis, angleEdges, histBinCentres);
                normDist = dirRel2Goal_histCounts_norm./totalDist_permute;
                sumNormDist = sum(normDist,1);
                normDistFactor = length(hdUp) ./ sumNormDist;
                normDistFactor = repmat(normDistFactor, length(histBinCentres), 1, 1);
                normDist = normDist .* normDistFactor;
                mrlNormShuffle = mrlRelDir(normDist, xAxis, yAxis, histBinCentres);
                mrlVal_Temp(s) = mrlNormShuffle.mrl;
            end
            
            %         mrlVal_raw = sort(mrlVal_raw, 'descend');
            %         mrlFocus(cellCounter).(trialTypes{tt}).raw.CI99P9 = mrlVal_raw(1);
            %         mrlFocus(cellCounter).(trialTypes{tt}).raw.CI95 = mrlVal_raw(50);
            %
            %         if mrlDataRaw.mrl > mrlVal_raw(50)
            %             signifCellsRaw.(trialTypes{tt}) = [signifCellsRaw.(trialTypes{tt}); ...
            %                 unitNames(c)];
            %         end
            
            mrlVal_norm = [mrlVal_norm; mrlVal_Temp];
            mrlVal_norm = sort(mrlVal_norm, 'descend');
            
            %%
%             if length(mrlVal_norm) == 1000
%                 htitle = [unitNames{c} '_' trialTypes{tt} '_relDirBootStrap'];
%                 h=figure('Name', htitle,'NumberTitle','off','Units',...
%                     'pixels', 'Position', [0 0 250 250], 'Color', [1 1 1], 'Visible', 'on');
%                 
%                 [counts, edges] = histcounts(mrlVal_norm);
%                 xAxis = edges(1:end-1) + (edges(2) - edges(1))/2;
%                 bar(xAxis, counts)
%                 hold on
%                 
%                 figCI95 = mrlVal_norm(50);
%                 figCI99P9 = mrlVal_norm(1);
%                 
%                 plot([figCI99P9, figCI99P9], [0, max(counts)], '-r')
%                 hold on
%                 
%                 plot([figCI95, figCI95], [0, max(counts)], '-r')
%                 hold on
%             
%                 plot([mrlDataNorm.mrl, mrlDataNorm.mrl], [0, max(counts)], '-k')
%                 hold on
%                 
%                 box off
%             
%                 set(gcf, 'PaperPositionMode', 'auto');
%                 set(gcf, 'InvertHardCopy', 'off')
%                 set(gcf,'color','w');
%                 set(gcf,'renderer','Painters')
%                 print('-depsc2', [htitle '.eps'])
%                 print('-djpeg', [htitle '.jpeg'])
%                 close
%             end
            %%
            
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
                         unitNames{c}];
                end
            end
        end
    end
end
%%
mrlFocus = reorderStructByUnit(mrlFocus);
% save mrlFocus_bursts_BONF.mat mrlFocus angleEdges histBinCentres % signifCellsNorm 
save mrlFocus_bursts_coarse.mat mrlFocus angleEdges histBinCentres signifCellsNorm 
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






