% compares original angle arrays with the corresponding best matching mean
% postures and calculates the deviation over time for each file in the
% specified directory

% load the representative postures matrix
% load(['/Users/abrown/Andre/wormVideos/robynTracking/optoG/zxIs6/' ...
%     'postures_50-centers_23-files_2000-framesPerFile.mat'])
load(['/Users/abrown/Andre/wormVideos/results-12-05-10/' ...
    'Laura Grundy/gene_NA/allele_NA/N2/on_food/XX/30m_wait/' ...
    'postures_150-centers_20-files_5000-framesPerFile.mat'])
postures = postures';

% set the root directory where the projected amplitudes are
% directory = ...
%     '/Users/abrown/Andre/wormVideos/robynTracking/optoG/zxIs6/'; % ljIs105 zxIs6
% directory = ...
%     '/Users/abrown/Andre/wormVideos/robynTracking/benzaldehyde/';
% directory = ['/Users/abrown/Andre/wormVideos/results-12-05-10/'...
%     'Laura Grundy/gene_NA/'];
directory = ['/Users/abrown/Andre/wormVideos/results-12-05-10/'...
    'Laura Grundy/unc-79/'];
% directory = ['/Users/abrown/Andre/wormVideos/results-12-05-10/' ...
%     'wild-isolates/'];

% re-create fileList looking for angle arrays
[fileList, ~] = dirSearch(directory, 'angleArray.mat');

% initialisation
devCell = cell(0);
r2Cell = cell(0);
r2Mat2 = NaN(numel(fileList), 1);
seqCell = cell(0);


% loop through the files in the directory
maxLength = 0;
for ii = 1:length(fileList)
    disp(ii/length(fileList))
    % load an angles file
    angles = cell2mat(struct2cell(load(fileList{ii}) ));
    
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
    
    % also get the R2 of the fit per worm using alternative function
    rSquared = calculateR2(anglesNoNaN, postures(stateSequence, :));
    
    % add deviation to cell
    devCell{ii} = D;
    r2Cell{ii} = r2;
    r2Mat2(ii) = rSquared;
    seqCell{ii} = stateSequence;
    if length(D) > maxLength
        maxLength = length(D);
    end
end

% get the mean deviation over the files and the corresponding standard
% deviation
devMat = NaN(length(fileList), maxLength);
for ii = 1:length(fileList)
    devMat(ii, 1:length(devCell{ii})) = devCell{ii};
end
seqMat = NaN(length(fileList), maxLength);
for ii = 1:length(fileList)
    seqMat(ii, 1:length(seqCell{ii})) = seqCell{ii};
end
r2Mat = NaN(length(fileList), maxLength);
for ii = 1:length(fileList)
    r2Mat(ii, 1:length(r2Cell{ii})) = r2Cell{ii};
end

meanDeviation = nanmean(devMat);
stdDeviation = nanstd(devMat);
semDeviation = nanstd(devMat) ./ sqrt(length(fileList));

% plot deviations
figure
patch([1:maxLength, maxLength:-1:1], ...
    [meanDeviation - semDeviation, ...
    meanDeviation(end:-1:1) + semDeviation(end:-1:1)], ...
    [0.7 0.7 0.7], 'EdgeColor', 'none')
line(1:maxLength, meanDeviation, 'LineWidth', 2, 'Color', [0.3 0.6 0.9])



meanR2 = nanmean(r2Mat);
stdR2 = nanstd(r2Mat);
semR2 = nanstd(r2Mat) ./ sqrt(length(fileList));

% plot deviations
figure
patch([1:maxLength, maxLength:-1:1], ...
    [meanR2 - stdR2, ...
    meanR2(end:-1:1) + stdR2(end:-1:1)], ...
    [0.7 0.7 0.7], 'EdgeColor', 'none')
line(1:maxLength, meanR2, 'LineWidth', 2, 'Color', [0.3 0.6 0.9])



% get the state densities over time
stateDensities = NaN(size(postures, 1), maxLength);
for ii = 1:maxLength
    stateDensities(:, ii) = histc(seqMat(:, ii), 1:size(postures, 1));
end

% plot the state densities
% imagesc(stateDensities)
clustergram(stateDensities * 0.8, ...
    'Cluster', 1, 'Linkage', 'complete', 'ColorMap', 'redbluecmap', ...
    'Symmetric', 'false')

