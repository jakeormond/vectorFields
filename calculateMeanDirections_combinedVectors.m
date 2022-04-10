% calculateMeanDirections_combinedVectorsV4
% FINAL ATTEMPT BEFORE RESUBMISSION TO GET GOAL 2 SINKS WHERE WE WANT THEM
% THIS IS THE ONE WE USE IN THE PAPER!!!!!!!!!!!!!!!!
%%
cd positionalData
load mazeCorners.mat
load goalCoordinates.mat
load platformLocations.mat
load platformOutlinesFinal.mat
load positionalDataByTrialType.mat
load directionMaps_PlatformsNEW.mat

cd ..

cd physiologyData
load burstSpikes.mat burstSpikes
load spikesByTrial.mat


cd direction
% load spikeHD_All_byPlatform.mat convergenceData
% load spikeHD_All_byPlatform_Shift_June2021V2.mat convergenceData
load mrlFocus_ctrlDistribution_coarse.mat
load mrlFocus_Shift_significantData.mat
% load mrlFocus_significantData_97P5.mat

trialTypes = fieldnames(pos);
%%
platOutlines(cellfun(@isempty,platOutlines)) = [];
platformPos = NaN(61, 2);
for p = 1:61
    platPos = cellfun(@(x) x{p}, platOutlines, 'un', 0);
    platPos = vertcat(platPos{:});
    platformPos(p,:) = mean(platPos,1);
end
%% get maze outline and goal coordinates in binned units for heat maps
binSize = 15;
xAxis = 1:binSize:frameSize(1);
yAxis = 1:binSize:frameSize(2);

mazeCornersX = interp1(xAxis, 1:length(xAxis), mazeCorners(:,1));
mazeCornersY = interp1(yAxis, 1:length(yAxis), mazeCorners(:,2));
mazeCorners4Maps = [mazeCornersX, mazeCornersY];
mazeCorners4Maps = [mazeCorners4Maps; mazeCorners4Maps(1,:)];

for tt = 1:length(trialTypes)
    goalCoorTemp = goalCoor.(['goal' num2str(goalID(tt))]);
    goalCoorX = interp1(xAxis, 1:length(xAxis), goalCoorTemp(:,1));
    goalCoorY = interp1(yAxis, 1:length(yAxis), goalCoorTemp(:,2));
    goalCoorTemp2 = [goalCoorX, goalCoorY];
    goalCoor4Maps{tt} = [goalCoorTemp2; goalCoorTemp2(1,:)];
    
    goalCoorX = interp1(xAxis, 1:length(xAxis), mean(goalCoorTemp(:,1)));
    goalCoorY = interp1(yAxis, 1:length(yAxis), mean(goalCoorTemp(:,2)));    
    goalBins{tt}(1) = round(goalCoorX);
    goalBins{tt}(2) = round(goalCoorY);
    
end

%%
dirEdges = 0:30:360;
radEdges = linspace(0, 2*pi, length(dirEdges));
radCentres = radEdges(1:end-1) + (radEdges(2) - radEdges(1))/2;

for tt = 1:length(trialTypes)
    cellCounter = 0;
    spikePlatsAll = [];
    for c = 1:length(mrlFocus)
        unit = mrlFocus(c).unit;
        units.(trialTypes{tt}){c,1} = unit;
        
        if isempty(mrlFocus(c).([trialTypes{tt} '_bursts'])) && ...
                isempty(mrlFocus(c).([trialTypes{tt} '_spikes']))
            continue
        end
        
        cellCounter = cellCounter + 1;
        
        spikePlats = [];
        spikeSamples = [];
        spikeHD = [];
        spikePos = [];
        
        for t = 1:length(pos.(trialTypes{tt}))            
            if ~isempty(mrlFocus(c).([trialTypes{tt} '_bursts']))
                type = 'bursts';
                mrlFocusTemp = mrlFocus(c).([trialTypes{tt} '_bursts']);
                
                if isempty(burstSpikes.(unit).(trialTypes{tt})(t).samples)
                    continue
                end                
                
                burstTemp = burstSpikes.(unit).(trialTypes{tt})(t);
                
%                 for b = 1:length(burstTemp.spikeHD)
%                     hdTemp = burstTemp.spikeHD{b};
%                     meanHD = rad2deg(circ_mean(deg2rad(hdTemp)));
%                     if meanHD < 0
%                         meanHD = meanHD + 360;
%                     elseif meanHD > 360
%                         meanHD = meanHD - 360;
%                     end
%                     burstTemp.spikeHD{b} = repmat(meanHD, length(hdTemp), 1); clear hdTemp
%                     
%                     posTemp = burstTemp.spikePos{b};
%                     meanPos = mean(posTemp,1);
%                     burstTemp.spikePos{b} = repmat(meanPos, size(posTemp,1), 1); clear posTemp
%                     
%                 end
                
                spikeSamplesTemp = vertcat(burstTemp.samples{:});
                spikeHDtemp = vertcat(burstTemp.spikeHD{:});
                spikePosTemp = vertcat(burstTemp.spikePos{:}); clear burstTemp
                                
%                 spikeSamplesTemp = vertcat(burstSpikes.(unit).(trialTypes{tt})(t).samples{:});
%                 spikeHDtemp = vertcat(burstSpikes.(unit).(trialTypes{tt})(t).spikeHD{:});
%                 spikePosTemp = vertcat(burstSpikes.(unit).(trialTypes{tt})(t).spikePos{:});
                
            else
                type = 'spikes';
                mrlFocusTemp = mrlFocus(c).([trialTypes{tt} '_spikes']);
                spInd = find(strcmp({spikesByTrial(:).unit}, unit));
                
                spikeSamplesTemp = spikesByTrial(spInd).(trialTypes{tt})(t).spikeSamples;
                spikeHDtemp = spikesByTrial(spInd).(trialTypes{tt})(t).spikeHD;
                spikePosTemp = spikesByTrial(spInd).(trialTypes{tt})(t).spikePos;
            end
                        
            platsTemp = vertcat(platformLocations.(trialTypes{tt}){t}.body);
            spikePlatsTemp = NaN(length(spikeHDtemp), 1);
            for s = 1:length(spikePlatsTemp)
                absDiff = abs(double(spikeSamplesTemp(s)) - pos.(trialTypes{tt})(t).sample);
                [~, ind2use] = min(absDiff);
                spikePlatsTemp(s) = platsTemp(ind2use);
            end
            
            spikePlats = [spikePlats; spikePlatsTemp];
            spikeSamples = [spikeSamples; spikeSamplesTemp];
            spikeHD = [spikeHD; spikeHDtemp];
            spikePos = [spikePos; spikePosTemp];
        end
        spikePlatsAll = [spikePlatsAll; spikePlats];
        %%
        dirMap = directionMaps.(trialTypes{tt}).headAngle;
        for p = 1:61
            platInd = find(spikePlats == p);
            if isempty(platInd)
                continue
            end
            spikeHD_byPlatform.(trialTypes{tt}).nSpikes(p,cellCounter) = ...
                length(platInd);
            
            platOcc = platOccupancy.(trialTypes{tt})(p);
            hdTemp = spikeHD(platInd);
            
            spikeHD_byPlatform.(trialTypes{tt}).meanRate(p,cellCounter) = ...
                length(hdTemp)/platOcc;
                
            meanDir = rad2deg(circ_mean(deg2rad(hdTemp)));
            if meanDir < 0
                meanDir = meanDir + 360;
            end
            
            spikeHD_byPlatform.(trialTypes{tt}).meanDir(p,cellCounter) = meanDir;
            spikeHD_byPlatform.(trialTypes{tt}).mrl(p,cellCounter) = ...
                circ_r(deg2rad(hdTemp));
            
            hdDist = histcounts(hdTemp, dirEdges);
            platHD = dirMap{p};
            platDist = histcounts(platHD, dirEdges);
            platTimeDist = platDist * (platOcc/sum(platDist));
            hdRatesDist = hdDist./platTimeDist;
            hdRatesDist(isnan(hdRatesDist)) = 0;
            hdRatesDist(isinf(hdRatesDist)) = 0;
            
            normMeanDir = rad2deg(circ_mean(radCentres', hdRatesDist'));
                            
            if normMeanDir < 0
                normMeanDir = normMeanDir + 360;
            end
            
            spikeHD_byPlatform.(trialTypes{tt}).meanDirNorm(p,cellCounter) = normMeanDir;
            spikeHD_byPlatform.(trialTypes{tt}).mrlNorm(p,cellCounter) = ...
                circ_r(radCentres', hdRatesDist');
            
            if isnan(spikeHD_byPlatform.(trialTypes{tt}).mrlNorm(p,cellCounter))
                jake = 1;
            end
        end
    end
    %% average vectors across cells
    vectors.(trialTypes{tt}).normDir = NaN(61,1);
    vectors.(trialTypes{tt}).normLength = NaN(61,1);
    
    nSpikesPerPlatform = sum(spikeHD_byPlatform.(trialTypes{tt}).nSpikes, 2);
       
    for p = 1:size(spikeHD_byPlatform.(trialTypes{tt}).meanRate, 1) 
        
        if nSpikesPerPlatform(p) < 100
            continue
        end
        
        normVectors = zeros(cellCounter, 2);
    
        for c = 1:cellCounter
            firingRate = spikeHD_byPlatform.(trialTypes{tt}).meanRate(p,c);
            mrlNorm = spikeHD_byPlatform.(trialTypes{tt}).mrlNorm(p,c);
            dirNorm = spikeHD_byPlatform.(trialTypes{tt}).meanDirNorm(p,c);
            
            scalFactorNorm = firingRate * mrlNorm;
            
            radDirNorm = deg2rad(dirNorm - 90);
            if radDirNorm < 0
                radDirNorm = radDirNorm + 2*pi;
            end            
            unitDist(1) = cos(radDirNorm);
            unitDist(2) = sin(radDirNorm);            
            normVectors(c,:) = unitDist * scalFactorNorm;
        end
        
        meanNormVector = nanmean(normVectors, 1);
        normVectorLength = sqrt(meanNormVector(1)^2 + meanNormVector(2)^2);
        if normVectorLength == 0
            continue
        end
        
        normDir = rad2deg(atan2(meanNormVector(2), meanNormVector(1))) + 90;
        if normDir < 0
            normDir = normDir + 360;
        end
        vectors.(trialTypes{tt}).normDir(p) = normDir;
        vectors.(trialTypes{tt}).normLength(p) = normVectorLength;
        
    end  
    
    %% calculate ConSink just from platform directions
    hdTemp = vertcat(vectors.(trialTypes{tt}).normDir(:));
    normLength = vertcat(vectors.(trialTypes{tt}).normLength(:));
    nSpikes = round(normLength*100);

    posTemp = platformPos;
    posTemp(isnan(hdTemp),:) = [];
    nSpikes(isnan(nSpikes)) = [];
    hdTemp(isnan(hdTemp)) = [];
    
    hdCell = {};
    posCell = {};
    for p = 1:length(hdTemp)
        hdCell{p,1} = repmat(hdTemp(p), nSpikes(p), 1);
        posCell{p,1} = repmat(posTemp(p,:), nSpikes(p), 1);    
    end
    
    hdTemp2 = vertcat(hdCell{:});
    posTemp2 = vertcat(posCell{:});
    mrlData = relativeDirectionFunction(posTemp2, hdTemp2, ... % posTemp and hdTemp for V5, 
        spikePlatsAll, relDirDists.(trialTypes{tt}), ... % posTemp2 and hdTemp2 for V6
        purDirDists.(trialTypes{tt}), angleEdges, frameSize);
    
    convData.(trialTypes{tt}) = mrlData.norm;
    
    goalMRL.(trialTypes{tt})(1) = ...
        mrlData.norm.allMRL(goalBins{1}(2), goalBins{1}(1));
%     goalMRL.(trialTypes{tt})(2) = ...
%         mrlData.norm.allMRL(goalBins{2}(2), goalBins{2}(1));
    
        
    
    %%
    if tt == 1
     
        htitle = ['meanDirections-Shift-vectors_June2021V6']; % v3 has averaged bursts, v4 does not
        h=figure('Name', htitle,'NumberTitle','off','Units',... % V5 has 1 direction per platform
            'pixels', 'Position', [50 50 750 500], 'Color', [1 1 1], ... % V6 has direction repeated on each platform
            'Visible', 'on');                                            % according to vector length       
    end
    
    subplot(2, 3, 3*(tt-1) + 1)
    
    mazeCornersTemp = [mazeCorners; mazeCorners(1, :)];
    plot(mazeCornersTemp(:,1), frameSize(2) - mazeCornersTemp(:,2), ...
        '-k', 'LineWidth', 2)
    
    hold on
    goalStr = ['goal' num2str(goalID(tt))];
    goalCoorTemp = [goalCoor.(goalStr); goalCoor.(goalStr)(1,:)];
    plot(goalCoorTemp(:,1), frameSize(2) - goalCoorTemp(:,2), '-k', 'LineWidth', 2)
    hold on
    
    % convPoint = convergenceData.(trialTypes{tt}).norm.coor;
    convPoint = mrlData.norm.coor;
    plot(convPoint(:,1), frameSize(2) - convPoint(:,2), ...
        '.r', 'MarkerSize', 20)
    hold on
    % dirCoor = convergenceData.(trialTypes{tt}).norm.dirCoor;
    dirCoor = mrlData.norm.dirCoor;
    plot(dirCoor(:,1), frameSize(2) - dirCoor(:,2), ...
        '.b', 'MarkerSize', 20)
    hold on
    
    for p = 1:61
        if isnan(vectors.(trialTypes{tt}).normDir(p))
            continue
        end
        
        position = platformPos(p,:);
        
        meanDir = vectors.(trialTypes{tt}).normDir(p);
        
        [u,v] = pol2cart(deg2rad(meanDir - 90), 35);
        
        xCoor = [position(1), position(1) - u];
        yCoor = [frameSize(2) - position(2), frameSize(2) - (position(2) + v)];
        
        plot(xCoor, yCoor, '-k', 'LineWidth', 1)
        hold on
        
        axis equal
        box off
    end
    %% direction map
    h1 = subplot(2, 3, 3*(tt-1) + 2);
    direction = rad2deg(convData.(trialTypes{tt}).allDir) + 180;
    direction = round(direction);
    direction(direction == 0) = 360;
    
    image(direction)
    cmap = flipud(hsv(360));
    % colormap(h1, hsv(360))
    colormap(h1, cmap)
    
    hold on
    plot(mazeCorners4Maps(:,1), mazeCorners4Maps(:,2), '-w', 'LineWidth', 2)
    
    hold on
    plot(goalCoor4Maps{tt}(:,1), goalCoor4Maps{tt}(:,2), '-w', 'LineWidth', 2)
    %% mrl map
    h2 = subplot(2, 3, 3*(tt-1) + 3);
    
    mrl = convData.(trialTypes{tt}).allMRL;
    mrl = (mrl./max(mrl(:)) * 63) + 1;
    
    image(mrl)
    colormap(h2, jet(64))
    
    hold on
    plot(mazeCorners4Maps(:,1), mazeCorners4Maps(:,2), '-w', 'LineWidth', 2)
    
    hold on
    plot(goalCoor4Maps{tt}(:,1), goalCoor4Maps{tt}(:,2), '-w', 'LineWidth', 2)
    
    box off
end


set(gcf, 'Renderer', 'Painters')
set(gcf, 'PaperPositionMode', 'auto');
print('-depsc2', [htitle '.eps'])
print('-djpeg', [htitle '.jpeg']);
% printeps(1, htitle)
close
%%
save spikeHD_Vectors_byPlatform_Shift_Jun2021V6.mat spikeHD_byPlatform ...
    vectors units convData goalMRL
    





