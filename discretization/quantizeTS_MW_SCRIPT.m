% go through entire data set to do some of the basic operations


% ----------------- parameters for run ------------------------------------
postureNum = 90;
warpNum = 1000;
minFrameNum = 100; % the minimum number of non-NaN frames allowed
flickThresh = 0;
verbose = 0;
ds = 5; % amount of downsampling, only every ds'th point is used
% -------------------------------------------------------------------------


% set the root directory where the projected amplitudes are
directory = '/Users/abrown/Andre/wormVideos/CeNDR/feature-files/';

% load postures
load([directory '_syntax-files/_CeNDR-postures_' num2str(postureNum) ...
    '-centers_198-files_1000-framesPerFile.mat'])
postures = postures';

% get the list of files (includes data from mutant classes)
[fileList, ~] = dirSearch(directory, 'angleArray.mat');

% drop data from videos with only 5 worms
dropInds = ~cellfun(@isempty, strfind(fileList, '_worms5_'));
fileList(dropInds) = [];

% get the worm names
strainNames = cell(numel(fileList), 1);
for ii = 1:numel(fileList)
    % get the filename
    [~, filename, ~] = fileparts(fileList{ii});
    
    % get the filename
    underPositions = strfind(filename, '_');
    
    % get the string with the strain name
    strainName = filename(1:underPositions(1)-1);
    
    % remove spaces and add to list
    strainNames{ii} = strrep(strainName, ' ', '');
end

% get the unique worm names
uniqueNames = unique(strainNames);


% loop through unique names
for kk = 1:numel(uniqueNames)
    % get the indices of the current strain
    currentNameInds = find(strcmp(strainNames, uniqueNames(kk)));
    
    % quantize each time series and save the results to disk
    for ii = 1:numel(currentNameInds)
        disp([kk/numel(uniqueNames), ii/numel(currentNameInds)])
        
%         try
            % load an angles file
            anglesTotal = cell2mat(struct2cell(load(fileList{currentNameInds(ii)}) ));
            
            % load the corresponding trajectories data
            [~, filename, ~] = fileparts(fileList{currentNameInds(ii)});
            trajFileName = [directory, strrep(filename, '_angleArray', '_featuresN.hdf5')];
            trajData = h5read(trajFileName, '/trajectories_data/');
            timeStamps = trajData.frame_number + 1;
            skelInds = trajData.skeleton_id + 1;
            wormIds = trajData.worm_index_joined + 1;
            uniqueWorms = unique(wormIds);
            
            % initialise
            seqCell = cell(0, 2);
            
            % loop through each tracked worm
            for jj = 1:numel(uniqueWorms)
                
                % get the indices for the current worm
                currentWormId = wormIds == uniqueWorms(jj);
                currentSkelInds = skelInds(currentWormId);
                
                % leave out zero indices which indicate dropped skeletons
                angles = ...
                    anglesTotal(currentSkelInds(currentSkelInds ~= 0), :);
                
                angles = angles(1:ds:end, :);
                
                % check number of good frames
                if sum(~isnan(angles(:, 1))) < minFrameNum/ds
                    continue
                end
                
                % drop any leading NaN values
                if isnan(angles(1, 1))
                    firstGoodPoint = find(~isnan(angles(:, 1)), 1, 'first');
                    angles(1:firstGoodPoint-1, :) = [];
                end
                
                % drop any trailing NaN values
                if isnan(angles(end, 1))
                    lastGoodPoint = find(~isnan(angles(:, 1)), 1, 'last');
                    angles(lastGoodPoint+1:end, :) = [];
                end
                
                % interpolate NaN gaps
                anglesNoNaN = angles;
                
                for mm = 1:size(angles, 2)
                    pAmp = angles(:, mm);
                    pAmp(isnan(pAmp)) = interp1(find(~isnan(pAmp)),...
                        pAmp(~isnan(pAmp)), find(isnan(pAmp)),'linear');
                    anglesNoNaN(:, mm) = pAmp;
                end
                
                % quantize the time series
                [stateSequence, ~] = ...
                    quantizeTS(anglesNoNaN, postures, warpNum, flickThresh);
                
                % add the current sequence to the sequence cell
                seqCell{jj, 1} = uniqueWorms(jj);
                seqCell{jj, 2} = stateSequence;
                
                % if in verbose mode, plot the overlap between the centers and the
                % original angle data for some examples
                if verbose
                    % re-calculate stateSequence without warping for plotting
                    stateSequence = quantizeTS(anglesNoNaN, postures, 1, flickThresh);
                    for mm = 1:length(stateSequence)
                        % plot the current skeleton angles
                        plot(1:size(anglesNoNaN, 2), anglesNoNaN(mm, :), '.', ...
                            'Color', [0.3 0.6 0.9])
                        hold on
                        % plot the best fitting cluster centre
                        line(1:size(postures, 2), postures(stateSequence(mm, 1), :), ...
                            'Color', [0.2 0.2 0.2], 'LineWidth', 2)
                        
                        ylim([-2.5 2.5])
                        hold off
                        
                        %     waitforbuttonpress
                        pause(0.1)
                        getframe;
                    end
                end
            end
            
            % drop any empty rows in seqCell
            dropInds = cellfun(@isempty, seqCell(:, 1));
            seqCell(dropInds, :) = [];
            
            % save the state sequence     
            seqFileName = strrep(fileList{currentNameInds(ii)}, '_angleArray', ...
                '_stateSequence_90-CeNDR-means_ds-5');
            save(seqFileName, 'seqCell')
           
%         catch
%             disp(['bad file name: ' fileList{currentNameInds(ii)}])
%         end
    end
end