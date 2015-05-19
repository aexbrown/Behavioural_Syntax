% go through entire data set to do some of the basic operations

postureNum = 90;

% load the representative postures matrix
% load(['/Users/abrown/Andre/wormVideos/results-12-05-10/' ...
%     'wild-isolates/postures_50-centers_200-files_500-framesPerFile.mat'])
load(['/Users/abrown/Andre/wormVideos/results-12-05-10/'...
    'Laura Grundy/gene_NA/allele_NA/N2/on_food/XX/30m_wait/' ...
    'postures_' num2str(postureNum) '-centers_20-files_5000-framesPerFile.mat'])
postures = postures';
distList = pdist(postures);

% set the root directory where the projected amplitudes are
directory = ...
    '/Users/abrown/Andre/wormVideos/robynTracking/';
% directory = ['/Users/abrown/Andre/wormVideos/results-12-05-10/'...
%     'Laura Grundy/gene_NA/allele_NA/N2/on_food/XX/30m_wait/'];
% directory = ['/Users/abrown/Andre/wormVideos/results-12-05-10/'...
%     'wild-isolates/CB4856/'];

% re-create fileList looking for angle arrays
[fileList, ~] = dirSearch(directory, 'angleArray.mat');

% initialise
maxDist = zeros(size(postures, 1), 1);
tightness = NaN(numel(fileList), 1);


% ----------------- parameters for run ------------------------------------
warpNum = 1000;
minFrameNum = 150; % the minimum number of non-NaN frames allowed
flickThresh = 0;
verbose = 0;
ds = 5; % amount of downsampling, only every ds'th point is used
% -------------------------------------------------------------------------

% quantize each time series and save the results to disk
for ii = 1:length(fileList)
    try
    % load an angles file
    angles = cell2mat(struct2cell(load(fileList{ii}) ));
    angles = angles(1:ds:end, :);
    
    % check number of good frames
    if sum(~isnan(angles(:, 1))) < minFrameNum/ds
        continue
    end
    
%     % if the worm is in the "left" folder, invert all the amplitudes
%     if ~isempty(strfind(fileList{ii}, 'CCW')) || ~isempty(strfind(fileList{ii}, 'ccw'))
%         angles = angles * -1;
%     end
    % if the worm is in the "left" folder, invert all the amplitudes
    if ~isempty(strfind(fileList{ii}, '/L/'))
        angles = angles * -1;
    end
    
    % drop any leading NaN values
    if isnan(angles(1, 1))
        firstGoodPoint = find(~isnan(angles(:, 1)), 1, 'first');
        angles(1:firstGoodPoint-1, :) = [];
    end
    
    % drop any trailing NaN values
    if isnan(angles(1, 1))
        lastGoodPoint = find(~isnan(angles(:, 1)), 1, 'last');
        angles(lastGoodPoint+1:end, :) = [];
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
    
    
    % quantize the time series
    [stateSequence, distVec] = ...
        quantizeTS(anglesNoNaN, postures, warpNum, flickThresh);
    
    % get the distance between each cluster centre and its most distant
    % cluster member
    for jj = 1:size(postures, 1)
        currentMax = max(distVec(stateSequence(:, 1) == jj));
        if currentMax > maxDist(jj)
            maxDist(jj) = currentMax;
        end
    end
    
    % see how many of the inter-cluster distances are larger than the
    % median maxDist
    tightness(ii) = sum(distList > median(maxDist))/length(distList);
    disp(tightness(ii))
    
    % if in verbose mode, plot the overlap between the centers and the
    % original angle data for some examples
    if verbose
        % re-calculate stateSequence without warping for plotting
        stateSequence = quantizeTS(anglesNoNaN, postures, 1, flickThresh);
        for jj = 1:10:length(stateSequence)
            % plot the current skeleton angles
            plot(1:size(anglesNoNaN, 2), anglesNoNaN(jj, :), '.', ...
                'Color', [0.3 0.6 0.9])
            hold on
            % plot the best fitting cluster centre
            line(1:size(postures, 2), postures(stateSequence(jj, 1), :), ...
                'Color', [0.2 0.2 0.2], 'LineWidth', 2)
            
            ylim([-2.5 2.5])
            hold off
            
            %     waitforbuttonpress
            pause(0.1)
            getframe;
        end
    end
    
    % save the state sequence
    statesFileName = ['stateSequence_' num2str(postureNum) 'means_N2-' num2str(warpNum) ...
        '_ds-' num2str(ds)];
    statesFileName = ...
        strrep(fileList{ii}, 'angleArray', statesFileName);
    save(statesFileName, 'stateSequence')
    
    disp([num2str(ii/numel(fileList)*100) '% complete.'])
    catch
        disp(['bad file number: ' num2str(ii)])
    end
end

% save the maxDist vector in the parent directory
save([directory 'maxDist_' num2str(postureNum) 'means_N2-' num2str(warpNum) ...
    '_ds-' num2str(ds) '.mat'], 'maxDist')