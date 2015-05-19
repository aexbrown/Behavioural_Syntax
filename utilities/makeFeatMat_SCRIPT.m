% use Ev's histogram function to extract the 726 feature measures for a set
% of feature files.  This saves histogram files to the same locations as
% the feature files, but these can be removed once the script is finished
% if you're only interested in the summary statistics.

% set directory
directory = ['/Users/abrown/Andre/wormVideos/results-12-05-10/' ...
    'wild-isolates/'];

% find feature files
[fileList, ~] = dirSearch(directory, 'features.mat');

% initialise feature matrix
featMat = NaN(numel(fileList), 726);
tic;
% loop through the files in fileList
for ii = 1:numel(fileList)
    disp(ii/numel(fileList))
    % get (and make) the current filenames
    filename = fileList{ii};
    statsFile = strrep(filename, 'features.mat', 'stats.mat');
    histFile = strrep(filename, 'features.mat', 'histogram.mat');
    
    % get the feature summary statistics
    [statStruct, statsInfo] = features2StatsArray(statsFile, histFile, filename);
    
    % add the feature summary vector to featMat
    featMat(ii, :) = [statStruct.mean];
end
toc;

% remove any rows or columns that are all NaN
nanCols = find(all(isnan(featMat)));
nanRows = find(all(isnan(featMat), 2));
featMat(:, nanCols) = [];
featMat(nanRows,:) = [];
nameList = fileList;
nameList(nanRows) = [];

% remove feature names for all-NaN features
featNames = {statsInfo.name};
featNames(nanCols) = [];

% hack the worm names
wormNames = cell(numel(fileList), 1);
for ii = 1:numel(fileList)
    % get the positions of forward slashes in file name
    slashPositions = strfind(fileList{ii}, '/');
    
    % get the strain name
    wormNames{ii} = fileList{ii}(slashPositions(7)+1:slashPositions(8)-1);
end

% impute NaN values to the column means
meanVals = nanmean(featMat);
featMatNoNaN = featMat;
for ii = 1:size(featMat, 1)
    % get the NaN values in the current row
    nanVals = isnan(featMat(ii, :));
    
    % replace the nanVals
    featMatNoNaN(ii, nanVals) = meanVals(nanVals);
end

% export the feature matrix as a csv file
csvwrite('featMat.csv', featMat)
csvwrite('featMatNoNaN.csv', featMatNoNaN)

% normalise a plot featMat
meanMatrix = repmat(nanmean(featMatNoNaN), size(featMatNoNaN, 1), 1);
stdMatrix = repmat(nanstd(featMatNoNaN), size(featMatNoNaN, 1), 1);
featMatNorm = (featMatNoNaN - meanMatrix) ./ stdMatrix;

% export the normalised feature matrix
csvwrite('featMatNorm.csv', featMatNorm)
