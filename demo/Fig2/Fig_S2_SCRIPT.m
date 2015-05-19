% script to plot n-gram frequency distributions with fits
%
% Reproduces Fig. S2

% should repeats be included?
repeats = false;

% set the root directory
directory = ['/Users/abrown/Andre/wormVideos/results-12-05-10/'...
    'Laura Grundy/gene_NA/allele_NA/N2/on_food/XX/30m_wait/'];

% get the file names
[fileList, ~] = ...
    dirSearch(directory, 'stateSequence_90means_N2-1000_ds-5.mat');

% get some random files
fileNum = 400;
rng(97987) % seed for reproducibility

% repeat for different random file combinations
nMax = 5;
numIterations = 10;

% initialisation
seqCell = cell(numIterations, fileNum);
freqCell = cell(nMax, numIterations);

% set the colour matrix
colVec = rand(nMax, 3);
colVec = colVec ./ (repmat(sum(colVec, 2), 1, 3) * 0.8);
colVec(colVec > 1) = 1;

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
    
    % plot frequency distributions
    for ii = 2:nMax
        % make total n-gram sets for the real data and simulated sequences
        nGramsTotal = [];
        for kk = 1:numel(randInds)
            % get the current sequence and a shuffled version
            currentSeq = seqCell{jj, kk};
            
            % calculate the n-grams
            % check that length isn't too short in each case
            if length(currentSeq) > 100
                nGramsTotal = ...
                    [nGramsTotal; n_gramsNumerical(currentSeq, ii)];
            end
        end
        
        % get the counts
        [~, counts] = countUniqueRows(nGramsTotal);
        
        % calculate the rank distribution and plot it
        freq = counts / size(nGramsTotal, 1);
        freqCell{ii, jj} = freq;
    end
end



% plot individual runs of data with the mean of the simulation runs
edges = 0:1e-5:1e-2;
for ii = 2:nMax
    figure
    for jj = 1:numIterations
        data = freqCell{ii, jj};
%         data = random('wbl', 1e-5, 0.5, 1, 1e4);
        counts = histc(data, edges);
        line(edges, counts / trapz(edges, counts), ...
            'Color', colVec(ii, :))
        
        % get the maximum likelihood estimate of the distribution
        % parameters for exponential and Weibull distributions
        paramsExp = fitdist(data', 'exp');
        paramsWbl = fitdist(data', 'wbl');
        
        % add lines to plot
        line(edges, pdf(paramsExp, edges), 'Color', 'r')
        line(edges, pdf(paramsWbl, edges), 'Color', [0 0 0])
    end
%     set(gca, 'XScale', 'log')
%     set(gca, 'YScale', 'log')
%     set(gca,'ticklength',2*get(gca,'ticklength'))
%     xlim([1e-5, 1e-2])
%     ylim([1, 1e5])
    box on
end
