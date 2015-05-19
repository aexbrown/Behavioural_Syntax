% script to go through feature files and save angle arrays calculated from
% the skeleton data

% set the directory to search for feature files
directory = '/Users/abrown/Andre/wormVideos/celine/';

% get the list of feature files
[fileList, ~] = dirSearch(directory, 'features.mat');

% loop through the file list
for ii = 1:numel(fileList)
    disp(ii/numel(fileList))
    
    % get (and make) the current filenames
    filename = fileList{ii};
    
    % get the directory name
    anglesFile = strrep(filename, 'features.mat', 'angleArray.mat');
    meansFile = strrep(filename, 'features.mat', 'meanAngles.mat');
%     slashPositions = strfind(filename, '/');
%     anglesFile = [filename(1:slashPositions(end)) 'angleArray.mat'];
%     meansFile = [filename(1:slashPositions(end)) 'meanAngles.mat'];
    
    featData = cell2mat(struct2cell(load( fileList{ii}, ...
        'worm'  )));
    
    % calculate the angle array from the skeleton data
    [angleArray, meanAngles] = ...
        makeAngleArray(featData.posture.skeleton.x', ...
        featData.posture.skeleton.y');
    
    % save the angle array
    save(anglesFile, 'angleArray')
    save(meansFile, 'meanAngles')
end