% script to plot posture probabilities over time in optogenetics
% experiments
% 
% Reproduces Fig. 3D/E

% load postures
load(['/Users/abrown/Andre/wormVideos/results-12-05-10/' ...
    'Laura Grundy/gene_NA/allele_NA/N2/on_food/XX/30m_wait/' ...
    'postures_90-centers_20-files_5000-framesPerFile.mat'])
postures = postures';

% set the root directory where the projected amplitudes are  ljIs105 zxIs6
wormName = 'ljIs105';
directory = ...
    ['/Users/abrown/Andre/wormVideos/robynTracking/optoG/' wormName '/'];

fps = 30; % the frame rate of the movies
ds = 5; % amount of downsampling, only every ds'th point is used
fps = fps/ds; % the frame rate after downsampling

% re-create fileList looking for angle arrays
[fileList, ~] = dirSearch(directory, 'angleArray.mat');

% initialisation
r2Cell = cell(0);
seqCell = cell(0);

% loop through the files in the directory
maxLength = 0;
minLength = Inf;
for ii = 1:length(fileList)
    disp(ii/length(fileList))
    % load an angles file
    angles = cell2mat(struct2cell(load(fileList{ii}) ));
    angles = angles(1:ds:end, :);
    
    % if the worm is in the "left" folder, invert all the amplitudes
    if ~isempty(strfind(fileList{ii}, '/L/'))
        angles = angles * -1;
    end
    
    % interpolate NaN gaps
    anglesNoNaN = angles;
    
    % interpolate over NaN values
    for jj = 1:size(angles, 2)
        pAmp = angles(:, jj);
        pAmp(isnan(pAmp)) = interp1(find(~isnan(pAmp)),...
            pAmp(~isnan(pAmp)), find(isnan(pAmp)),'linear');
        anglesNoNaN(:, jj) = pAmp;
    end
    
    % find the distance between each posture and its nearest neighbour in
    % the matrix of representative postures
    [stateSequence, D] = knnsearch(anglesNoNaN, postures, 1);
    
    % convert D to R2 per frame using deviation of each posture (note, mean
    % angle of each post
    D_tot = sum(( anglesNoNaN - ...
        repmat(nanmean(anglesNoNaN, 2), 1, size(anglesNoNaN, 2)) ).^2, 2);
    
    r2 = 1 - (D ./ D_tot);
    
    % add deviation to cell
    r2Cell{ii} = r2;
    seqCell{ii} = stateSequence;
    if length(D) > maxLength
        maxLength = length(D);
    end
    if length(D) < minLength
        minLength = length(D);
    end
end

% get the state sequence and r2 values over time
seqMat = NaN(length(fileList), maxLength);
for ii = 1:length(fileList)
    seqMat(ii, 1:length(seqCell{ii})) = seqCell{ii};
end
r2Mat = NaN(length(fileList), maxLength);
for ii = 1:length(fileList)
    r2Mat(ii, 1:length(r2Cell{ii})) = r2Cell{ii};
end


% plot the fit quality over time
meanR2 = nanmean(r2Mat);
stdR2 = nanstd(r2Mat);
semR2 = nanstd(r2Mat) ./ sqrt(length(fileList));

% plot R2 over time
figure
patch([1:maxLength, maxLength:-1:1], ...
    [meanR2 - stdR2, ...
    meanR2(end:-1:1) + stdR2(end:-1:1)], ...
    [0.7 0.7 0.7], 'EdgeColor', 'none')
line(1:maxLength, meanR2, 'LineWidth', 2, 'Color', [0.3 0.6 0.9])
xlim([20, 45] * fps)
ylim([0.6, 1])

% get the state densities over time
stateDensities = NaN(size(postures, 1), maxLength);
for ii = 1:maxLength
    stateDensities(:, ii) = ...
        histc(seqMat(:, ii), 1:size(postures, 1)) ./ size(seqMat, 1);
end

% smooth the state densities
filtSize = 9;
sdSmooth = NaN(size(stateDensities));
for ii = 1:90
    sdSmooth(ii, :) = ...
        conv(stateDensities(ii, :), ones(filtSize, 1)/filtSize, 'same');
end

% plot the state densities.  Cluster to aid visualisation.
clg = clustergram(sdSmooth, 'Cluster', 1, 'Linkage', 'complete');

% it's easier to control the colour using image than clustergram, so take
% row labels and re-plot
rowOrder = str2num(char(clg.RowLabels));

% plot the entire sequence
figure
imagesc([0, size(sdSmooth, 2)] / fps, [1, 90], ...
    sdSmooth(rowOrder(end:-1:1), :), [0, 0.06])
cMap1 = [1:-0.1:0; 1:-0.1:0; ones(1, 11)]';
cMap2 = colormap('jet');
cMap = [cMap1; cMap2(9:end, :)];
colormap(cMap)

% make another figure and plot only the first stimulus portion
figure
imagesc([0, size(sdSmooth, 2)] / fps, [1, 90], ...
    sdSmooth(rowOrder(end:-1:1), :), [0, 0.06])
colormap(cMap)
xlim([20, 45])

% also plot a version that is averaged across the three repeated stimuli
sdAverage = ( sdSmooth(:, round(20*fps):round(45*fps)) + ...
    sdSmooth(:, round(55*fps):round(80*fps)) + ...
    sdSmooth(:, round(90*fps):round(115*fps)) ) / 3;
figure
imagesc([0, size(sdAverage, 2)] / fps, [1, 90], ...
    sdAverage(rowOrder(end:-1:1), :), [0, 0.06]) 
colormap(cMap)




% repeat clustering and plotting using a row-normalised matrix
sdAverageRowMean = mean(sdAverage, 2);
sdAverageRowStd = std(sdAverage, [], 2);
sdAverageNorm = ...
    (sdAverage - repmat(sdAverageRowMean, 1, size(sdAverage, 2))) ./ ...
    repmat(sdAverageRowStd, 1, size(sdAverage, 2));

sdAverageNorm(isnan(sdAverageNorm(:, 1)), :) = 0;

% plot the state densities.  Cluster to aid visualisation.
clg = clustergram(sdAverageNorm, 'Cluster', 1, 'Linkage', 'complete');

% it's easier to control the colour using image than clustergram, so take
% row labels and re-plot
rowOrder = str2num(char(clg.RowLabels));
save(['rowOrder_' wormName '.mat'], 'rowOrder')

figure
imagesc([0, size(sdAverage, 2)] / fps, [1, 90], ...
    sdAverageNorm(rowOrder(end:-1:1), :), [-2.5, 2.5]) 
cmap = cbrewer('div', 'RdBu', 50); % RdYlBu
colormap(cmap(end:-1:1, :))


