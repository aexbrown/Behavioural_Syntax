% find a set of representative postures by loading angleArrays from a
% directory and performing k-means clustering.

% clear angleTotal so that it doesn't accumulate from previous script runs
% in case workspace was not cleared
clear angleTotal

% set the directory to search for feature files
directory = '/Users/abrown/Andre/wormVideos/CeNDR/feature-files/_syntax-files/';

% ---------------------------- Set Parameters -----------------------------
minFrames = 10000; % the minimum number of acceptable frames in a file
% fileNum = 20; !!! fileNum is set below in the loop !!!
frameNum = 1000; % the number of random frames to take per file
k = [5 10:10:200 225 250 300 400 500]; % the number of centres to use for k-means clustering
% -------------------------------------------------------------------------

% seed random number generator for reproducibility
rng(0954);

% get the list of files (includes data from mutant classes)
[fileList, ~] = dirSearch(directory, 'angleArray.mat');

% drop data from videos with only 5 worms
dropInds = ~cellfun(@isempty, strfind(fileList, '_worms5_'));
fileList(dropInds) = [];

% hack the worm names
strainNames = cell(numel(fileList), 1);
for ii = 1:numel(fileList)
    % get the filename
    [~, filename, ~] = fileparts(fileList{ii});
    
    % get the filename
    underPositions = strfind(filename, '_');
    
    % get the string with the strain or mutant/allele name
    strainNames{ii} = filename(1:underPositions(1)-1);
end

% get the unique worm names
uniqueNames = unique(strainNames);

angleTotal = [];

% loop through unique names
for jj = 1:numel(uniqueNames)
    disp(['Angle collection ' num2str(jj/numel(uniqueNames)*100) '% complete.'])
    
    % reset fileNum.  N.B. if fileNum is >= to the number of files in
    % fileList, then the script will simply run deterministically on all
    % files in the list
    fileNum = 1;
    
    % get the indices of the current strain
    currentNameInds = find(strcmp(strainNames, uniqueNames(jj)));
    
    isRandom = 1;
    % check if random files or all files should be analysed
    if fileNum >= numel(currentNameInds)
        fileNum = numel(currentNameInds);
        isRandom = 0;
    end
    
    % loop through files
    for ii = 1:fileNum
        if isRandom
            % get a random file
            % N.B. if directory has a small number of files there could be
            % repeats keep searching until angleArray has enough elements
            angleArray = [];
            while size(angleArray, 1) < minFrames || sum(~isnan(angleArray(:, 1))) < 2
                fileInd = ...
                    currentNameInds(randi(numel(currentNameInds), 1));
                
                % load the angle array data
                angleArray = cell2mat(struct2cell(load(fileList{fileInd})));
            end
        else
            fileInd = currentNameInds(ii);
            angleArray = cell2mat(struct2cell(load(fileList{fileInd})));
            if sum(~isnan(angleArray(:, 1))) < 2
                continue
            end
        end
        
        % take a random subset of frames (could include repeats). Take 20%
        % more than requested in case of NaNs
        frameInds = randi(size(angleArray, 1), round(frameNum * 1.2), 1);
        
        % exclude any selected NaN rows
        nanRows = find(isnan(angleArray(:, 1)));
        frameInds = setdiff(frameInds, nanRows);
        
        % get only requested number of frames
        if numel(frameInds) < frameNum
            disp(['Warning: only ' num2str(numel(frameInds)) ...
                ' frames included. File contained more than 20% NaN frames.'])
        else
            frameInds = frameInds(1:frameNum);
        end
        
        % add the current angle array to the total
        angleTotal = [angleTotal; angleArray(frameInds, :)];
    end
end

for ii = 1:numel(k)
    disp(['Clustering ' num2str(ii/numel(k)*100) '% complete.'])
    
    % do the clustering
    [~, postures] = kmeans(angleTotal', k(ii));
    
    % save the representative postures in a mat file
    save([directory '_CeNDR-postures_' num2str(k(ii)) '-centers_' ...
        num2str(fileNum * numel(uniqueNames)) ...
        '-files_' num2str(frameNum) '-framesPerFile.mat'], 'postures')
end
