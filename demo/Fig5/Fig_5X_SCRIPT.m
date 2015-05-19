% script to check whether N2 postures also fit wild-isolates
%
% Reproduces Fig. 5X

% load the representative postures matrix
load(['/Users/abrown/Andre/wormVideos/results-12-05-10/' ...
    'Laura Grundy/gene_NA/allele_NA/N2/on_food/XX/30m_wait/' ...
    'postures_90-centers_20-files_5000-framesPerFile.mat'])
postures = postures';

% set amount of downsampling, only every ds'th point is used
ds = 5;

% set the root directory where the projected amplitudes are
directory = ['/Users/abrown/Andre/wormVideos/results-12-05-10/' ...
    'wild-isolates/'];
% get the file list
[fileList, ~] = dirSearch(directory, 'angleArray.mat');

% remove AQ2947 which is actually just the CGC N2 from this analysis
dropInds = zeros(numel(fileList), 1);
for ii = 1:numel(fileList)
    if ~isempty(strfind(fileList{ii}, '/AQ2947/'))
        dropInds(ii) = 1;
    end
end
fileList(logical(dropInds)) = [];

% hack the worm names
wormNames = cell(numel(fileList), 1);
for ii = 1:numel(fileList)
    % get the positions of forward slashes in file name
    slashPositions = strfind(fileList{ii}, '/');
    
    % get the strain name
    wormNames{ii} = fileList{ii}(slashPositions(7)+1:slashPositions(8)-1);
end

% get the unique worm names
uniqueNames = unique(wormNames);

% initialisation
r2Cell = cell(0);
r2Hist = zeros(201, 1);
seqCell = cell(0);


% loop through the files in the directory
maxLength = 0;
for jj = 1:numel(uniqueNames)
    disp(jj/numel(uniqueNames))
    
    % get the current file indices
    currentInds = find(strcmp(wormNames, uniqueNames{jj}));
    
    for ii = 1:numel(currentInds)
        % load an angles file
        angles = cell2mat(struct2cell(load(fileList{currentInds(ii)}) ));
        angles = angles(1:ds:end, :);
        
        % skip files that are all NaN
        if sum(~isnan(angles(:))) == 0
            continue
        end
        
        % if the worm is in the "left" folder, invert all the amplitudes
        if ~isempty(strfind(fileList{ii}, '/L/'))
            angles = angles * -1;
        end
        
        % interpolate NaN gaps
        anglesNoNaN = angles;
        
        % interpolate over NaN values
        for kk = 1:size(angles, 2)
            pAmp = angles(:, kk);
            pAmp(isnan(pAmp)) = interp1(find(~isnan(pAmp)),...
                pAmp(~isnan(pAmp)), find(isnan(pAmp)),'linear');
            anglesNoNaN(:, kk) = pAmp;
        end
        
        % find the distance between each posture and its nearest neighbour
        % in the matrix of representative postures
        [stateSequence, D] = knnsearch(anglesNoNaN, postures, 1);
        
        % convert D to R2 per frame using deviation of each posture (note,
        % mean angle of each post
        D_tot = sum(( anglesNoNaN - ...
            repmat(nanmean(anglesNoNaN, 2), 1, size(anglesNoNaN, 2))...
            ).^2, 2);
        
        r2 = 1 - (D ./ D_tot);
        
        % get distribution of r2 values
        counts = histc(r2, 0:0.005:1);
        r2Hist = r2Hist + counts;
        
        % add deviation to cell
        r2Cell{ii} = r2;
        seqCell{ii} = stateSequence;
        if length(D) > maxLength
            maxLength = length(D);
        end
    end
    % add a line for the current histogram
    line(0:0.005:1, r2Hist/sum(r2Hist), ...
        'LineWidth', 1, 'Color', [0.7 0.7 0.7])
end
