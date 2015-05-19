% script to check whether on-food postures can also fit data in other
% conditions (e.g. off food or during chemotaxis)
% 
% Reproduces Fig. 4B

% load the representative postures matrix
postureNum = 90;
load(['/Users/abrown/Andre/wormVideos/results-12-05-10/' ...
    'Laura Grundy/gene_NA/allele_NA/N2/on_food/XX/30m_wait/' ...
    'postures_' num2str(postureNum) '-centers_20-files_5000-framesPerFile.mat'])
postures = postures';

% set the root directory where the projected amplitudes are
directory = ...
    '/Users/abrown/Andre/wormVideos/robynTracking/benzaldehyde/'; %benzaldehyde no-food
% directory = ['/Users/abrown/Andre/wormVideos/results-12-05-10/'...
%     'Laura Grundy/gene_NA/allele_NA/N2/on_food/XX/30m_wait/'];
% directory = ['/Users/abrown/Andre/wormVideos/results-12-05-10/'...
%     'Laura Grundy/unc-79/'];
% directory = ['/Users/abrown/Andre/wormVideos/results-12-05-10/' ...
%     'wild-isolates/'];

ds = 5; % amount of downsampling, only every ds'th point is used
edges = 0:0.05:1; % edges to use in plotting histogram

% get the file list
[fileList, ~] = dirSearch(directory, 'angleArray.mat');

% for N2 data, get a random subset of 100 files
if numel(fileList) > 500
    rng(10)
    randInds = randi(numel(fileList), 100, 1);
    fileList = fileList(randInds);
end

% initialisation
r2Cell = cell(postureNum, 1);
r2Hist = zeros(size(edges, 2), 1);


% loop through the files in the directory
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
    
    % convert D to R2 per frame using deviation of each posture
    D_tot = sum(( anglesNoNaN - ...
        repmat(nanmean(anglesNoNaN, 2), 1, size(anglesNoNaN, 2)) ).^2, 2);
    
    r2 = 1 - (D ./ D_tot);
    
    % get distribution of r2 values
    counts = histc(r2, edges);
    r2Hist = r2Hist + counts;
    
    % record the coefficient of determination separately for each posture
    for kk = 1:postureNum
        % get the indices of the current posture
        currentInds = stateSequence == kk;
        
        r2Cell{kk} = [r2Cell{kk}; r2(currentInds)];
    end
end

% plot the histogram of r2 values
line(edges, r2Hist/sum(r2Hist), 'LineWidth', 2, 'Color', [0.3 0.9 0.1])


% plot the mean and standard deviation for each posture
colorVec = [0.3 0.9 0.1];
figure
hold on
for ii = 1:postureNum
%     errorbar(ii+0.25, mean(r2Cell{ii}), std(r2Cell{ii}), 'Color', [0.1 0.3 0.9])
%     plot(ii+0.25, quantile(r2Cell{ii}, 0.15), 'o', 'Color', 'r')
    line([ii, ii]-0.5, [mean(r2Cell{ii}) - std(r2Cell{ii}), ...
        mean(r2Cell{ii}) + std(r2Cell{ii})], ...
        'LineWidth', 2, 'Color', colorVec)
    plot(ii-0.5, mean(r2Cell{ii}), ...
        'o', 'Color', colorVec)
end
hold off

