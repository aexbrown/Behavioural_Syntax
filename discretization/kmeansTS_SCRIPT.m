% find a set of representative postures by loading angleArrays from a
% directory and performing k-means clustering.

% clear angleTotal so that it doesn't accumulate from previous script runs
% in case workspace was not totally cleared
clear angleTotal

% set the root directory
% directory = ['/Users/abrown/Andre/wormVideos/results-12-05-10/'...
%     'Laura Grundy/'];
% directory = ['/Users/abrown/Andre/wormVideos/results-12-05-10/'...
%     'Laura Grundy/gene_NA/allele_NA/N2/on_food/XX/30m_wait/'];
% directory = ['/Users/abrown/Andre/wormVideos/results-12-05-10/' ...
%     'wild-isolates/'];
directory = ...
    '/Users/abrown/Andre/wormVideos/robynTracking/benzaldehyde/';

% ---------------------------- Set Parameters -----------------------------
minFrames = 1000; % the minimum number of acceptable frames in a file
fileNum = 200; % the number of random files to process.  N.B. if fileNum is
               % >= to the number of files in fileList, then the script
               % will simply run deterministically on all files in the list
frameNum = 5000; % the number of random frames to take per file
k = 500; % the number of centres to use for k-means clustering
% -------------------------------------------------------------------------

% seed random number generator for reproducibility
rng(49230);

% get the list of files (includes data from mutant classes)
[fileList, ~] = dirSearch(directory, 'angleArray.mat');

isRandom = 1;
% check if random files or all files should be analysed
if fileNum >= numel(fileList)
    fileNum = numel(fileList);
    isRandom = 0;
end

% loop through files
angleTotal = [];
for ii = 1:fileNum
    if isRandom
        % get a random file
        % N.B. if directory has a small number of files there could be
        % repeats keep searching until angleArray has enough elements
        angleArray = [];
        while size(angleArray, 1) < minFrames || sum(~isnan(angleArray(:, 1))) < 2
            fileInd = randi(numel(fileList), 1);
            
            % load the angle array data
            angleArray = cell2mat(struct2cell(load(fileList{fileInd})));
        end
    else
        fileInd = ii;
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
    
    % since extrapolation was not used, there could still be some remaining
    % NaNs.  Remove these here.
    nanRows = find(isnan(anglesNoNaN(:, 1)));
    anglesNoNaN(nanRows, :) = [];
    
    % take a random subset of frames (could include repeats)
    frameInds = randi(size(anglesNoNaN, 1), frameNum, 1);
    anglesNoNaN = anglesNoNaN(frameInds, :);
    
    % if the worm is in the "left" folder, invert all the angles.
    if ~isempty(strfind(fileList{fileInd}, '/L/'))
        anglesNoNaN = anglesNoNaN * -1;
    end
    
    % add the current angle array to the total
    angleTotal = [angleTotal; anglesNoNaN];
end

% do the clustering
[~, postures] = kmeans(angleTotal', k);

% save the representative postures in a mat file
save([directory 'BENZALDEHYDE_postures_' num2str(k) '-centers_' num2str(fileNum) ...
    '-files_' num2str(frameNum) '-framesPerFile.mat'], 'postures')
