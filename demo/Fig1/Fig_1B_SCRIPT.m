% script to plot posture distribution
% 
% Reproduces Fig. 1B

% set the root directory
directory = ['/Users/abrown/Andre/wormVideos/results-12-05-10/'...
    'Laura Grundy/gene_NA/allele_NA/N2/on_food/XX/30m_wait/'];

% get the file names
[fileList, ~] = ...
    dirSearch(directory, 'stateSequence_90means_N2-1000_ds-5.mat');

% get some random files
fileNum = 20;
rng(53) % seed for reproducibility

randInds = randi(numel(fileList), fileNum, 1);

% loop through files to make larger and larger strings
stateSequence = [];
for ii = 1:numel(randInds)
    % import the data
    stateData = cell2mat(struct2cell(load(fileList{randInds(ii)})));
    stateSequenceCurrent = stateData(:, 1)';
    stateSequence = [stateSequence stateSequenceCurrent];
end
    
% get the distribution of postures, sorted from most common
[stateCount, ~] = histc(stateSequence, 1:90);
[~, sortInds] = sort(stateCount, 'descend');

% plot the state probabilities
figure
bar(stateCount(sortInds) ./ sum(stateCount), 'BarStyle', 'hist')
shading flat

% load the postures
load(['/Users/abrown/Andre/wormVideos/results-12-05-10/' ...
    'Laura Grundy/gene_NA/allele_NA/N2/on_food/XX/30m_wait/' ...
    'postures_90-centers_20-files_5000-framesPerFile.mat'])
postures = postures';

% select postures that span the range of frequencies
selectedPostures = sortInds([1, 16, 26, 35, 58, 69, 85, 90]);

% plot some sample postures
for ii = 1:numel(selectedPostures)
    % reconstruct the skeleton coordinates
    [x, y] = angle2skel(postures(selectedPostures(ii), :)', 0, 1);
    
    % plot the skeleton and add a point for the head
    figure
    plot(x, y, 'LineWidth', 12, 'Color', 'b')
    axis equal
    xlim([-0.25, 1.25])
    hold on
    plot(x(1), y(1), '.', 'MarkerSize', 35, 'Color', 'r')
    hold off
end




% also make a supplementary figure which shows all of the postures
% plot some sample postures
figure
[xoffset, yoffset] = meshgrid(1:9, 1:10);
for ii = 1:numel(sortInds)
    % reconstruct the skeleton coordinates
    [x, y] = angle2skel(postures(sortInds(ii), :)', 0, 1);
    
    % plot the skeleton and add a point for the head
    plot(x + 1.5*xoffset(ii), y - 0.5*yoffset(ii), ...
        'LineWidth', 4, 'Color', 'b')
    hold on
    plot(x(1) + 1.5*xoffset(ii), y(1) - 0.5*yoffset(ii), ...
        '.', 'MarkerSize', 10, 'Color', 'r')
end
xlim([1, 15])
axis equal
