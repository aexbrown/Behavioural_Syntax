% script to plot n-gram frequency distributions and fit to stretched
% exponentials
%
% Reproduces Fig. 5A (Zipf plot)

% should repeats be included?
repeats = false;

% set the root directory
directory = ['/Users/abrown/Andre/wormVideos/results-12-05-10/'...
    'wild-isolates/'];

% get the file names
[fileList, ~] = ...
    dirSearch(directory, 'stateSequence_90means_N2-1000_ds-5.mat');

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

% size of n-gram to use
n = 3;

% initialisation and model set up
fitMat = NaN(numel(uniqueNames), 3); % fit results for each strain
modelfun = @(b,x)(b(1)*exp( -(x/b(2)).^b(3) )); % the stretched exponential
dropNum = 0; % number of points to neglect from start of freqCurve for fit
betaStart = [0.03, 1, 0.2]; % starting parameters for the nonlinear fit
seqCell = cell(numel(uniqueNames), 1); % cell for storing state sequences
freqCell = cell(numel(uniqueNames), 1); % frequency distributions for each strain
pct99 = NaN(numel(uniqueNames), 1); % the 99th percentile for each curve
rank99 = NaN(numel(uniqueNames), 1); % the rank corresponding to each pct99
top1Frac = NaN(numel(uniqueNames), 1); % the fraction of time spent using the top 1% of the repertoire


% set the colour matrix
colVec = rand(numel(uniqueNames), 3);
colVec = colVec ./ (repmat(sum(colVec, 2), 1, 3) * 0.8);
colVec(colVec > 1) = 1;

for jj = 1:numel(uniqueNames)
    %     disp(jj/numel(uniqueNames))
    
    % get the file names of the current strain
    currentInds = find(strcmp(wormNames, uniqueNames{jj}));
    
    % loop through files to make larger and larger strings
    stateSequence = [];
    for ii = 1:numel(currentInds)
        % import the data
        stateData = cell2mat(struct2cell(load(fileList{currentInds(ii)})));
        % ds is 1 for expansion because although data are downsampled, the
        % data used to train the language model were also downsampled
        if repeats
            expandedSeq = expandSequence(stateData, 1);
        else
            expandedSeq = stateData;
        end
        stateSequence = [stateSequence expandedSeq(:, 1)'];
    end
    seqCell{jj} = stateSequence;
    
    % get n-grams
    nGrams = n_gramsNumerical(stateSequence, n);
    
    % get the counts
    [~, counts] = countUniqueRows(nGrams);
    
    % sort the n-grams by frequency
    [~, sortInds] = sort(counts, 'descend');
    
    % calculate the rank distribution and plot it
    freqCurve = counts(sortInds) / size(nGrams, 1);
    freqCell{jj} = freqCurve;
    
    % fit the frequency distribution
    beta = nlinfit(1:(length(freqCurve)-dropNum), ...
        freqCurve(1 + dropNum:end), modelfun, betaStart);
    fitMat(jj, :) = beta;
    
    % get the 99th percentile of the frequency distribution and the
    % corresponding rank
    pct99(jj) = quantile(freqCurve, 0.99);
    rank99(jj) = ceil(length(freqCurve)*0.01);
    top1Frac(jj) = sum(freqCurve(1:rank99(jj)));
    
    % plot frequency curve
%     figure
    line(1:length(freqCurve), freqCurve, 'Color', colVec(jj, :))
%     % add the fit with mean parameters over the top
%     a = modelfun(beta, 1:length(freqCurve));
%     line(1:length(freqCurve), a, 'LineWidth', 2, 'Color', [0.3, 0.6, 0.9])
%     text(100, 0.01, ['tau = ' num2str(beta(2)/beta(3) * gamma(1/beta(3)))])
end

% adjust plot scale/ticks
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
ylim([2e-6, 4e-2])
set(gca,'ticklength',2*get(gca,'ticklength'))

% % add the fit with mean parameters over the top
% a = modelfun(median(fitMat), 1:length(freqCurve));
% line(1:length(freqCurve), a, 'LineWidth', 2, 'Color', [0.3, 0.6, 0.9])
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')


% plot lines for the mean ± std of rank99
line([mean(rank99), mean(rank99)], [1e-6, 1e-1], 'Color', 'r')
line([mean(rank99), mean(rank99)] - std(rank99), [1e-6, 1e-1], 'Color', 'r')
line([mean(rank99), mean(rank99)] + std(rank99), [1e-6, 1e-1], 'Color', 'r')







