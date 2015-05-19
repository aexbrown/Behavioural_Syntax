% interpolate NaN values from feature files

% set directory
% directory = '/Users/abrown/Andre/wormVideos/robynTracking/';
directory = ['/Users/abrown/Andre/wormVideos/results-12-05-10/part/'...
    'gene_NA/allele_NA/N2/on_food/XX/30m_wait/'];

% find feature files
[fileList, ~] = dirSearch(directory, 'features.mat');

% interpolate linearly over NaN values, but drop leading and trailing NaN
% values
trim = 1;
interpNaN(fileList, 'projectedAmpsNoNaN.mat', 100, trim, 1)
