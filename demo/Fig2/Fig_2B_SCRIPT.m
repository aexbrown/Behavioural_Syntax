% script to plot n-gram frequency distributions
%
% Reproduces Fig. 2B (Zipf plot)

% should repeats be included?
repeats = false;

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

% get some random files
fileNum = 400;
fileNumSim = 400;
rng(53) % seed for reproducibility

% repeat for different random file combinations
nMax = 5;
numIterations = 10;

% initialisation
seqCell = cell(numIterations, fileNum); % cell for storing state sequences
seqCellSim = cell(numIterations, fileNum); % cell for storing state sequences
freqCell = cell(nMax - 1, 1);
freqCellShuf = cell(nMax - 1, 1);
freqCellSim = cell(nMax - 1, 1);
pct99 = NaN(nMax, numIterations, 1); % the 99th percentile for each curve
rank99 = NaN(nMax, numIterations, 1); % the rank corresponding to each pct99
top1Frac = NaN(nMax, numIterations, 1); % the fraction of time spent using the top 1% of the repertoire


% set the colour matrix
colVec = rand(nMax, 3);
colVec = colVec ./ (repmat(sum(colVec, 2), 1, 3) * 0.8);
colVec(colVec > 1) = 1;

minSizeVec = Inf(1, nMax - 1);
for jj = 1:numIterations
    disp(jj/numIterations)
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
    
    % also get random subset of simulated sequences without repeats
    permInds2 = randperm(size(simSeqMat, 1));
    randInds2 = permInds2(1:fileNumSim);
    for ii = 1:numel(randInds2)
        % get a row from simSeqMat
        simSeqPart = simSeqMat(randInds2(ii), :);
        
        % drop zeros introduced by dlmread for padding
        simSeqPart = simSeqPart(simSeqPart ~= 0);
        seqCellSim{jj, ii} = simSeqPart;
    end
    
    % make Zipf plots for different n n-grams
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
            % check that length isn't too short in each case
            if length(currentSeq) > 100
                nGramsTotal = ...
                    [nGramsTotal; n_gramsNumerical(currentSeq, ii)];
                nGramsTotalShuf = ...
                    [nGramsTotalShuf; n_gramsNumerical(currentSeqShuf, ii)];
            end
            
            if length(seqCellSim{jj, kk}) > 100
                nGramsTotalSim = ...
                    [nGramsTotalSim; n_gramsNumerical(seqCellSim{jj, kk}, ii)];
            end
        end
        
        % get the counts
        [~, counts] = countUniqueRows(nGramsTotal);
        [~, countsShuf] = countUniqueRows(nGramsTotalShuf);
        [~, countsSim] = countUniqueRows(nGramsTotalSim);
        
        % calculate the rank distribution and plot it
        freqCurve = sort(counts / size(nGramsTotal, 1), 'descend');
        freqCurveShuf = sort(countsShuf / size(nGramsTotalShuf, 1), 'descend');
        freqCurveSim = sort(countsSim / size(nGramsTotalSim, 1), 'descend');
        
        % get the 99th percentile of the frequency distribution and the
        % corresponding rank
        pct99(ii, jj) = quantile(freqCurve, 0.99);
        rank99(ii, jj) = ceil(length(freqCurve)*0.01);
        top1Frac(ii, jj) = sum(freqCurve(1:rank99(ii, jj)));
        
        % for subsequent average, we only data up to the shortest
        % accumulation curve, so keep track of min length
        if length(freqCurve) < minSizeVec(ii - 1)
            minSizeVec(ii - 1) = length(freqCurve);
        end
        if length(freqCurveSim) < minSizeVec(ii - 1)
            minSizeVec(ii - 1) = length(freqCurveSim);
        end
        
        % also add curves to cells for processing outsite loop.  On first
        % iteration add, on subsequent iterations concatenate
        if jj == 1
            freqCell{ii - 1} = freqCurve;
            freqCellShuf{ii - 1} = freqCurveShuf;
            freqCellSim{ii - 1} = freqCurveSim;
        else
            freqCell{ii - 1} = ...
                [freqCell{ii - 1}(:, 1:minSizeVec(ii - 1)); ...
                freqCurve(1:minSizeVec(ii - 1))];
            freqCellShuf{ii - 1} = ...
                [freqCellShuf{ii - 1}(:, 1:minSizeVec(ii - 1)); ...
                freqCurveShuf(1:minSizeVec(ii - 1))];
            freqCellSim{ii - 1} = ...
                [freqCellSim{ii - 1}(:, 1:minSizeVec(ii - 1)); ...
                freqCurveSim(1:minSizeVec(ii - 1))];
        end
    end
end


% plot individual runs of data with the mean of the simulation runs
figure
for ii = 2:nMax
    for jj = 1:numIterations
        line(1:minSizeVec(ii-1), freqCell{ii - 1}(jj, :), ...
            'Color', colVec(ii, :))
    end
    meanCurveSim = mean(freqCellSim{ii - 1});
    line(1:length(meanCurveSim), meanCurveSim, ...
        'LineWidth', 2, 'Color', [0 0 0])
end
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca,'ticklength',2*get(gca,'ticklength'))
xlim([1, 1e5])
ylim([1e-6, 2e-2])

% plot lines for the mean ± std of rank99
line([mean(rank99(3, :)), mean(rank99(3, :))], ...
    [1e-6, 1e-1], 'Color', 'r')
line([mean(rank99(3, :)), mean(rank99(3, :))] - std(rank99(3, :)), ...
    [1e-6, 1e-1], 'Color', 'r')
line([mean(rank99(3, :)), mean(rank99(3, :))] + std(rank99(3, :)), ...
    [1e-6, 1e-1], 'Color', 'r')




