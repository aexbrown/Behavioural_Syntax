% go through entire data set to do some of the basic operations


% ----------------- parameters for run ------------------------------------
postureNum = 90;
warpNum = 1000;
minFrameNum = 10000; % the minimum number of non-NaN frames allowed
flickThresh = 0;
verbose = 0;
ds = 5; % amount of downsampling, only every ds'th point is used
% -------------------------------------------------------------------------


% set the root directory where the projected amplitudes are
% directory = ...
%     '/Users/abrown/Andre/wormVideos/robynTracking/';
directory = '/Users/abrown/Andre/wormVideos/results-12-05-10/';
% directory = ['/Users/abrown/Andre/wormVideos/results-12-05-10/'...
%     'Laura Grundy/npr-1/'];
% directory = ['/Users/abrown/Andre/wormVideos/results-12-05-10/'...
%     'wild-isolates/CB4856/'];
postureDir = '/Users/abrown/Andre/wormVideos/results-12-05-10/postures/';


% get the list of files (includes data from mutant classes)
[fileList, ~] = dirSearch(directory, 'angleArray.mat');
fileStructPostures = dir(postureDir);

% exclude "no_wait" data
dropInds = [];
for ii = 1:numel(fileList)
    if ~isempty(strfind(fileList{ii}, 'no_wait'))
        dropInds = [dropInds, ii];
    end
end
fileList(dropInds) = [];

% hack the worm names
wormNames = cell(numel(fileList), 1);
for ii = 1:numel(fileList)
    % check for N2 hermaphrodites
    if ~isempty(strfind(fileList{ii}, '/N2/'))
        if ~isempty(strfind(fileList{ii}, '/XO/'))
            wormNames{ii} = 'N2_male';
        else
            wormNames{ii} = 'N2_herm';
        end
        continue
    end
    
    % get the positions of forward slashes in file name
    slashPositions = strfind(fileList{ii}, '/');
    
    % get the string with the strain or mutant/allele name
    wormName = fileList{ii}(slashPositions(7)+1:slashPositions(9)-1);
    
    %  remove '/on_food/' (present in wild isolate names and replace
    % slashes with underscores
    wormName = strrep(wormName, '/on_food', '');
    wormName = strrep(wormName, '/', '_');
    wormNames{ii} = wormName;
end

% get the unique worm names
uniqueNames = unique(wormNames);


% loop through unique names
for kk = 1:numel(uniqueNames)    
    for ii = 1:numel(fileStructPostures)
        % check for a match with the current strain name
        if ~isempty(strfind(fileStructPostures(ii).name, uniqueNames{kk})) && ...
            ~isempty(strfind(fileStructPostures(ii).name, ...
                ['postures_' num2str(postureNum)]))
            % load the postures file
            load([postureDir, fileStructPostures(ii).name])
            postures = postures';
            continue
        end
    end
    
    % get the indices of the current strain
    currentNameInds = find(strcmp(wormNames, uniqueNames(kk)));
    
    % quantize each time series and save the results to disk
    for ii = 1:numel(currentNameInds)
        disp([kk/numel(uniqueNames), ii/numel(currentNameInds)])

        try
            % load an angles file
            angles = cell2mat(struct2cell(load(fileList{currentNameInds(ii)}) ));
            angles = angles(1:ds:end, :);
            
            % check number of good frames
            if sum(~isnan(angles(:, 1))) < minFrameNum/ds
                continue
            end
            
            %     % if the worm is in the counter-clockwise folder, invert all the amplitudes
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
            
            for jj = 1:size(angles, 2)
                pAmp = angles(:, jj);
                pAmp(isnan(pAmp)) = interp1(find(~isnan(pAmp)),...
                    pAmp(~isnan(pAmp)), find(isnan(pAmp)),'linear');
                anglesNoNaN(:, jj) = pAmp;
            end
            
            
            % quantize the time series
            [stateSequence, distVec] = ...
                quantizeTS(anglesNoNaN, postures, warpNum, flickThresh);
            
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
            statesFileName = [num2str(currentNameInds(ii)) ...
                '_stateSequence_' num2str(postureNum) 'means_' ...
                uniqueNames{kk} '-' num2str(warpNum) ...
                '_ds-' num2str(ds)];
            statesFileName = ...
                strrep(fileList{currentNameInds(ii)}, 'angleArray', statesFileName);
            save(statesFileName, 'stateSequence')
        catch
            disp(['bad file name: ' fileList{currentNameInds(ii)}])
        end
    end
end