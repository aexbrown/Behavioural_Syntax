% script to go through feature files and save angle arrays calculated from
% the skeleton data.
% 
% Adapted to work with Tierpsy Tracker data

% set the directory to search for feature files
directory = '/Users/abrown/Andre/wormVideos/CeNDR/feature-files/';

% get the list of feature files
[fileList, ~] = dirSearch(directory, 'featuresN.hdf5');

% loop through the file list
for ii = 1:numel(fileList)
    disp(ii/numel(fileList))
    
    % read features file with skeletons
    skeletons = h5read(fileList{ii}, '/coordinates/skeletons');
    
    % calculate the angle array from the skeleton data
    [angleArray, meanAngles] = ...
        makeAngleArray(squeeze(skeletons(1, :, :))', ...
        squeeze(skeletons(2, :, :))');
    
    % make output filename
    [~, filename, ~] = fileparts(fileList{ii});
    anglesFile = strrep(filename, 'featuresN', 'angleArray.mat');
    meansFile = strrep(filename, 'featuresN', 'meanAngles.mat');
    
    % save the angle array
    save([directory '_syntax-files/' anglesFile], 'angleArray')
    save([directory '_syntax-files/' meansFile], 'meanAngles')
end