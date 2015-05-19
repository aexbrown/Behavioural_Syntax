% script to plot skeletons, their best fitting template, the state
% sequence, and duration of each state.
%
% Reproduces Fig. 1C

% load the features, angles, and templates
% the chosen worm is in an 'R' directory, so no need to flip it
load('/Users/abrown/Andre/wormVideos/results-12-05-10/Laura Grundy/gene_NA/allele_NA/N2/on_food/XX/30m_wait/R/tracker_5/2011-05-12___12_46_32/N2 on food L_2011_05_12__12_46_32__7_features.mat')
load('/Users/abrown/Andre/wormVideos/results-12-05-10/Laura Grundy/gene_NA/allele_NA/N2/on_food/XX/30m_wait/R/tracker_5/2011-05-12___12_46_32/angleArray.mat')
load('/Users/abrown/Andre/wormVideos/results-12-05-10/Laura Grundy/gene_NA/allele_NA/N2/on_food/XX/30m_wait/R/tracker_5/2011-05-12___12_46_32/meanAngles.mat')
load('/Users/abrown/Andre/wormVideos/results-12-05-10/Laura Grundy/gene_NA/allele_NA/N2/on_food/XX/30m_wait/postures_90-centers_20-files_5000-framesPerFile.mat')
% transpose postures
postures = postures';

% choose frames to plot
frames = (130:5:189)*1;

% choose downsampling rate
ds = 5;
angleArray = angleArray(1:ds:end, :);
meanAngles = meanAngles(1:ds:end);
skelXArray = worm.posture.skeleton.x(:, 1:ds:end);
skelYArray = worm.posture.skeleton.y(:, 1:ds:end);

% interpolate NaN gaps
anglesNoNaN = angleArray;

% interpolate over NaN values
for jj = 1:size(angleArray, 2)
    pAmp = angleArray(:, jj);
    pAmp(isnan(pAmp)) = interp1(find(~isnan(pAmp)),...
        pAmp(~isnan(pAmp)), find(isnan(pAmp)),'linear');
    anglesNoNaN(:, jj) = pAmp;
end

% get the state sequence and deviation between originals and templates
[stateSequence, ~] = knnsearch(anglesNoNaN, postures, 1);





% plot the original skeletons and their nearest neighbours
offsetFactor = 100;
figure
hold on
% loop through the frames
for ii = 1:length(frames)
    frame = frames(ii);
    
    % get the skeleton
    skelX = skelXArray(:, frame);
    skelY = skelYArray(:, frame);
    
    % get the next non-NaN frame
    while isnan(skelX(1))
        frame = frame + 1;
        
        % get the skeleton
        skelX = skelXArray(:, frame);
        skelY = skelYArray(:, frame);
        
        % if this is the last frame of the sequence, we need to extend for
        % plotting below
        if ii == length(frames)
            frames(end) = frames(end) + 1;
        end
    end
    
    
    % offset the skeleton
    skelX = skelX - mean(skelX) + (frame - frames(1)) * offsetFactor;
    skelY = skelY - mean(skelY);
    arclength = sum(sqrt(diff(skelX).^2 + diff(skelY).^2));
    
    % plot the original skeleton
    plot(skelX, skelY, '.', 'MarkerSize', 20, 'Color', [0 0 0])
    
    % convert the current matching posture to xy-coordinates
    [matchX, matchY] = angle2skel(postures(stateSequence(frame), :)', ...
        meanAngles(frame), arclength);
    
    % plot the matching posture over the original skeleton
    plot(matchX - mean(matchX) + (frame - frames(1)) * offsetFactor, ...
        matchY - mean(matchY), 'Color', 'r', 'LineWidth', 3)
    
    % add a vertical line to indicate where the skeletons are in time
    plot([mean(skelX), mean(skelX)], [-600, -1000], 'Color', [0 0 0])
    
    % make x and y scales the same
    axis equal
end

% add the state sequence as text with bars indicating durations
subSequence = stateSequence(frames(1):frames(end));
jumpInds = diff([subSequence; -1]) ~= 0;
durations = diff([0; find(jumpInds)]);

statesNoRepeats = subSequence(jumpInds);

breakPoints = [0; cumsum(durations)*offsetFactor];
switchColor = 1;
for ii = 1:numel(breakPoints)-1
    switchColor = switchColor * -1;
    % add a line to indicate state duration
    if switchColor == 1
        line([breakPoints(ii), breakPoints(ii+1)] - offsetFactor/2, ...
            [-1000, -1000], 'LineWidth', 1, 'Color', [0.9 0.6 0.3])
    else
        line([breakPoints(ii), breakPoints(ii+1)] - offsetFactor/2, ...
            [-1000, -1000], 'LineWidth', 1, 'Color', [0.3 0.6 0.9])
    end
    % write state sequence as text
    text((breakPoints(ii) + breakPoints(ii+1))/2, ...
        -1200, num2str(statesNoRepeats(ii)), 'HorizontalAlignment', 'center')
end
hold off

