% script to plot n-gram accumulation curves 
% 
% Reproduces Fig. 2A

% should repeats be included?
repeats = false;

% what step size should be used for calculating the accumulation curve?
stepSize = 5000;

% set the root directory
directory = ['/Users/abrown/Andre/wormVideos/results-12-05-10/'...
    'Laura Grundy/gene_NA/allele_NA/N2/on_food/XX/30m_wait/'];

% get the file names
[fileList, ~] = ...
    dirSearch(directory, 'stateSequence_90means_N2-1000_ds-5.mat');
fileNameSim = ['/Users/abrown/Andre/wormVideos/roland/andre_20150310/' ...
    'lm_N2_onfood_90_no-repeats_k3/random_seqs_k3.txt'];

% load simulated state sequence
simSeqMat = dlmread(fileNameSim);


% set number of random files to get
fileNum = 100;
fileNumSim = 175;
rng(542) % seed for reproducibility

% repeat for different random file combinations
nMax = 5;
numIterations = 20;

% intialise
accCurveCell = cell(nMax - 1, 1);
accCurveSimCell = cell(nMax - 1, 1);
seqCell = cell(numIterations, fileNum); % cell for storing state sequences
seqCellSim = cell(numIterations, fileNum); % cell for storing state sequences

% set the colour matrix
colVec = rand(nMax, 3);
colVec = colVec ./ (repmat(sum(colVec, 2), 1, 3) * 0.8);
colVec(colVec > 1) = 1;

minSize = Inf;
for jj = 1:numIterations
    disp(jj)
    % get random indices, but without repeats
    permInds = randperm(numel(fileList));
    randInds = permInds(1:fileNum);
    
    % loop through files to make larger and larger strings
    for ii = 1:numel(randInds)
        % import the data
        stateData = cell2mat(struct2cell(load(fileList{randInds(ii)})));
        % ds is 1 for expansion because although data are downsampled, the
        % data used to train the language model were also downsampled
        if repeats
            expandedSeq = expandSequence(stateData, 1);
        else
            expandedSeq = stateData;
        end
        seqCell{jj, ii} = expandedSeq(:, 1)';
    end
    
    % also get random subset of simulated sequences
    permInds2 = randperm(size(simSeqMat, 1));
    randInds2 = permInds2(1:fileNumSim);
    for ii = 1:numel(randInds2)
        % get a row from simSeqMat
        simSeqPart = simSeqMat(randInds2(ii), :);
        
        % drop zeros introduced by dlmread for padding
        simSeqPart = simSeqPart(simSeqPart ~= 0);
        seqCellSim{jj, ii} = simSeqPart;
    end
    
    % plot accumulation curves for different n n-grams
    for ii = 2:nMax        
        % make total n-gram sets for the real data and simulated sequences
        nGramsTotal = [];
        nGramsTotalShuf = [];
        nGramsTotalSim = [];
        for kk = 1:numel(randInds)
            % get the current sequence and a shuffled version
            currentSeq = seqCell{jj, kk};
            currentSeqShuf = currentSeq(randperm(numel(currentSeq)));
            
            % calculate the n-grams
            nGramsTotal = ...
                [nGramsTotal; n_gramsNumerical(currentSeq, ii)];
            nGramsTotalShuf = ...
                [nGramsTotalShuf; n_gramsNumerical(currentSeqShuf, ii)];
            % some of the simulated sequences are very short.  Make sure
            % there are at least n states in the sequence
            if size(seqCellSim{jj, kk}, 2) >= ii
                nGramsTotalSim = ...
                    [nGramsTotalSim; n_gramsNumerical(seqCellSim{jj, kk}, ii)];
            end
        end
        
        % calculate accumulation curves
        accCurve = nGramAccumulateStep(nGramsTotal, stepSize);
        accCurveShuf = nGramAccumulateStep(nGramsTotalShuf, stepSize);
        accCurveSim = nGramAccumulateStep(nGramsTotalSim, stepSize);
        
        % for subsequent average, we only data up to the shortest
        % accumulation curve, so keep track of min length
        if length(accCurve) < minSize
            minSize = length(accCurve);
        end
        if length(accCurveSim) < minSize
            minSize = length(accCurveSim);
        end
        
        % also add curves to cells for processing outsite loop.  On first
        % iteration add, on subsequent iterations concatenate
        if jj == 1
            accCurveCell{ii - 1} = accCurve;
            accCurveShufCell{ii - 1} = accCurveShuf;
            accCurveSimCell{ii - 1} = accCurveSim;
        else
            accCurveCell{ii - 1} = ...
                [accCurveCell{ii - 1}(:, 1:minSize); accCurve(1:minSize)];
            accCurveShufCell{ii - 1} = ...
                [accCurveShufCell{ii - 1}(:, 1:minSize); ...
                accCurveShuf(1:minSize)];
            accCurveSimCell{ii - 1} = ...
                [accCurveSimCell{ii - 1}(:, 1:minSize); ...
                accCurveSim(1:minSize)];
        end
    end
end

% plot individual runs of data with the mean of the simulation runs
figure
for ii = 2:nMax
    for jj = 1:numIterations
        line([1, (1:minSize-1)*stepSize], accCurveCell{ii - 1}(jj, 1:minSize), ...
            'Color', colVec(ii, :))
    end
    meanCurveShuf = mean(accCurveShufCell{ii - 1}, 1);
    line([1, (1:minSize-1)*stepSize], meanCurveShuf, ...
        'LineWidth', 2, 'Color', 'b')
    
    meanCurveSim = mean(accCurveSimCell{ii - 1}, 1);
    line([1, (1:minSize-1)*stepSize], meanCurveSim, ...
        'LineWidth', 2, 'Color', 'r')
end
xlim([0, 10.8e4])
ylim([0, 6.5e4])
set(gca,'LineWidth',1,'TickLength',[0.015 0.015]);