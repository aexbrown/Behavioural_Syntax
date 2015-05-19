% script to plot posture probabilities over time in optogenetics
% experiments
% 
% Reproduces Fig. 3C

% load postures
load(['/Users/abrown/Andre/wormVideos/results-12-05-10/' ...
    'Laura Grundy/gene_NA/allele_NA/N2/on_food/XX/30m_wait/' ...
    'postures_90-centers_20-files_5000-framesPerFile.mat'])
postures = postures';

% set the root directory where the projected amplitudes are  ljIs105 zxIs6
directory = ...
    '/Users/abrown/Andre/wormVideos/robynTracking/optoG/zxIs6/';

fps = 30; % the frame rate of the movies
ds = 5; % amount of downsampling, only every ds'th point is used
fps = fps/ds; % the frame rate after downsampling
n = 3; % the value of n for n-grams

% divide sequences into regions of stimulation and regions outside
% stimulation.  Use 5 seconds after stimulation ceases to account for
% response time
onInds = [30, 65, 100, 130] * fps;
offInds = [0, 40, 75, 110] * fps;

% get list of state sequences
[fileList, ~] = dirSearch(directory, ...
    ['stateSequence_90means_N2-1000_ds-' num2str(ds) '.mat']);

% initialisation
seqCell = cell(numel(fileList), 1);
uniqueNGrams = [];

% loop through the files in the directory
for ii = 1:numel(fileList)
    % load an angles file
    stateData = cell2mat(struct2cell(load(fileList{ii}) ));
    
    % expand the sequence to undo the warping
    expandedSeq = expandSequence(stateData, 1);
    seqCell{ii} = expandedSeq(:, 1)';
    
    % get n-grams
    nGrams = n_gramsNumerical(expandedSeq(:, 1)', n);
    
    % get the counts
    [uniqueRows, ~] = countUniqueRows(nGrams);
    
    % append current unique rows to total uniqueNGrams
    currentNGrams = [uniqueNGrams; uniqueRows];
    
    % re-calculate unique n-grams from larger set
    [uniqueNGrams, ~] = countUniqueRows(currentNGrams);
end

% compare n-gram densities outside of stimulation and during stimulation
nonStimFreq = zeros(numel(fileList), size(uniqueNGrams, 1));
stimFreq = zeros(numel(fileList), size(uniqueNGrams, 1));


for ii = 1:numel(fileList)
    % loop through intervals
    for jj = 1:numel(onInds)
        % get the sequences for the current non-stimulus intervals
        currentNonStimSeq = seqCell{ii}(offInds(jj)+1:onInds(jj)-1);
        
        % add counts of any observed n-grams
        for kk = 1:size(uniqueNGrams, 1)
            count = ...
                numel( strfind(currentNonStimSeq, uniqueNGrams(kk, :)) );
            nonStimFreq(ii, kk) = nonStimFreq(ii, kk) + ...
                count/size(currentNonStimSeq, 2);
        end
        
        % repeat for the stimulus intervals, unless it is the last interval
        % (the last interval is always non-stimulus, so there is no
        % corresponding stimulus interval)
        if jj < numel(onInds)
            currentStimSeq = seqCell{ii}(onInds(jj):offInds(jj+1));
            
            % add counts of any observed n-grams
            for kk = 1:size(uniqueNGrams, 1)
                count = ...
                    numel( strfind(currentStimSeq, uniqueNGrams(kk, :)) );
                stimFreq(ii, kk) = stimFreq(ii, kk) + ...
                    count/size(currentStimSeq, 2);
            end
        end
    end
end


% find any n-grams that are significantly differently used in one of the
% conditions using a ranksum test.  Record the effect size as well.
pVals = NaN(size(uniqueNGrams, 1), 1);
meanDiff = NaN(size(uniqueNGrams, 1), 1);
for kk = 1:size(uniqueNGrams, 1)
    % test for significance
    p = ranksum(stimFreq(:, kk), nonStimFreq(:, kk));
    pVals(kk) = p;
    
    % get the difference between the mean values
    meanDiff(kk) = mean(stimFreq(:, kk)) - mean(nonStimFreq(:, kk));
end

% estimate q-values
bhFDR = mafdr(pVals, 'BHFDR', true);


figure
hist(meanDiff(bhFDR < 0.05), 20)
xlim([-0.15, 0.15])

% find the n-grams that are significantly enriched during/following
% stimulation
upInds = find(bhFDR < 0.05 & meanDiff > 0);

% sort upInds based on FDR
[~, sortInds] = sort(bhFDR(upInds));
upIndsSorted = upInds(sortInds);

% plot the sequences that are most different
for ii = 1:numel(upIndsSorted)
    if ii == 10
        break
    end
    % plot the skeleton and add a point for the head
    figure
    plotSequence(postures, uniqueNGrams(upIndsSorted(ii), :))
end



% find the n-grams that are significantly enriched during/following
% stimulation
downInds = find(bhFDR < 0.05 & meanDiff < 0);

% sort upInds based on FDR
[~, sortInds] = sort(bhFDR(downInds));
downIndsSorted = downInds(sortInds);

% plot the sequences that are most different
for ii = 1:numel(downIndsSorted)
    if ii == 10
        break
    end
    % plot the skeleton and add a point for the head
    figure
    plotSequence(postures, uniqueNGrams(downIndsSorted(ii), :), 'b', 'r')
end

