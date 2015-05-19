% load feature matrices from a directory, perform mRMR feature selection
% and save the vectors of selected feature indices

% how many features should be selected?
numFeatures = 1000;

% set the directory and get the list of files
directory = ['/Users/abrown/Andre/wormVideos/results-12-05-10/' ...
    'analysisResults/'];
[fileList, ~] = dirSearch(directory, 'featMat_gTerminals*.mat');

% load the worm names
wormNames = struct2cell(load([directory, 'wormNames.mat']));
wormNames = wormNames{:};

% convert worm names to integers indicating grouped mutant classes
[~, ~, classInts] = unique(wormNames);

% loop through files
for ii = 2:numel(fileList)
    % get current file
    featMat = cell2mat(struct2cell(load(fileList{ii})));
    
    % select features
    tic;
    featureIndices = mRMR_feature_select(featMat, classInts, numFeatures);
    toc;
    
    % save selected features
    save(strrep(fileList{ii}, 'featMat', 'featInds'), 'featureIndices')
end