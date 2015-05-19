% get the average R2 for different total template number (k)
%
% Reproduces Fig. 1A

% set the root directory
directory = ['/Users/abrown/Andre/wormVideos/results-12-05-10/'...
    'Laura Grundy/gene_NA/allele_NA/N2/on_food/XX/30m_wait/'];

% ---------------------------- Set Parameters -----------------------------
minFrames = 15000; % the minimum number of acceptable frames in a file
fileNum = 50; % the number of random files to process.
frameNum = 5000; % the number of random frames to take per file
kVec = [5 10:10:200 225 250 300 400 500]; % the number of postures to use
% -------------------------------------------------------------------------

% initialise
R2cell = cell(length(kVec), 1);
seqcell = cell(length(kVec), 1);

% seed random number generator for reproducibility
rng(49230);

% get the list of files (includes data from mutant classes)
[fileList, ~] = dirSearch(directory, 'angleArray.mat');

% loop through posture numbers
for nn = 1:length(kVec)
    disp(nn/length(kVec))
    k = kVec(nn);
    
    % load the k template postures
    load([directory 'N2_postures_' ...
        num2str(k) '-centers_100-files_5000-framesPerFile.mat'])
    postures = postures';

    DTotal = [];
    R2Total = [];
    seqTotal = [];
    % loop through files
    for ii = 1:fileNum
        % get a random file
        % N.B. if directory has a small number of files there could be
        % repeats. Keep searching until angleArray has enough elements
        angleArray = 1;
        while sum(~isnan(angleArray(:, 1))) < minFrames
            fileInd = randi(numel(fileList), 1);
            
            % load the angle array data
            angleArray = cell2mat(struct2cell(load(fileList{fileInd})));
        end
        
        % instead of dropping NaN values, interpolate over NaNs
        % initialise anglesNoNaN
        anglesNoNaN = angleArray;
        
        % interpolate over NaN values
        for jj = 1:size(angleArray, 2)
            pAmp = angleArray(:, jj);
            pAmp(isnan(pAmp)) = interp1(find(~isnan(pAmp)),...
                pAmp(~isnan(pAmp)), find(isnan(pAmp)),'linear');
            anglesNoNaN(:, jj) = pAmp;
        end
        
        % since extrapolation was not used, there could still be some
        % remaining NaNs.  Remove these here.
        nanRows = find(isnan(anglesNoNaN(:, 1)));
        anglesNoNaN(nanRows, :) = [];
        
        % take a random subset of frames (could include repeats)
        frameInds = randi(size(anglesNoNaN, 1), frameNum, 1);
        anglesNoNaN = anglesNoNaN(frameInds, :);
        
        % if the worm is in the "left" folder, invert all the angles.
        if ~isempty(strfind(fileList{fileInd}, '/L/'))
            anglesNoNaN = anglesNoNaN * -1;
        end
        
        % get the nearest neighbours of the current angle array
        [stateSequence, ~] = knnsearch(anglesNoNaN, postures, 1);
        seqTotal = [seqTotal, stateSequence];
        
        % get the correlation coefficient
        R2Vec = NaN(length(stateSequence), 1);
        for jj = 1:length(stateSequence)
            R = corrcoef(anglesNoNaN(jj, :), ...
                postures(stateSequence(jj), :));
            R2Vec(jj) = R(1, 2)^2;
        end
        R2Total = [R2Total; R2Vec];
    end
    
    % add R2 for this posture number to the output cell
    R2cell{nn} = R2Total;
    seqcell{nn} = seqTotal;
end

% get the mean R2 for each file and value of k
meanMat = [];
for ii = 1:length(kVec)
    % get the mean of each set
    meanVec = [];
    for jj = 1:frameNum:length(R2cell{ii})
        meanVec = [meanVec; mean(R2cell{ii}(jj:jj + frameNum - 1))];
    end
    meanMat = [meanMat meanVec];
end

% plot the mean ± std or sem for the R2 values
meanR2 = mean(meanMat);
stdR2 = std(meanMat);

% plot deviations
figure
patch([kVec, kVec(end:-1:1)], ...
    [meanR2 - stdR2, ...
    meanR2(end:-1:1) + stdR2(end:-1:1)], ...
    [0.7 0.7 0.7], 'EdgeColor', 'none')
line(kVec, meanR2, 'LineWidth', 2, 'Color', [0.3 0.6 0.9])

figure
errorbar(kVec, meanR2, stdR2, '.', ...
    'MarkerSize', 13, 'Color', [0.3 0.6 0.9])


% also make a plot of the rate of mean R2 change
plot(kVec(2:end), diff(meanR2)./diff(kVec), '-')
a = diff(meanR2)./diff(kVec);
line([0, 500], [mean(a(9:end)), mean(a(9:end))])
xlim([0, 520])