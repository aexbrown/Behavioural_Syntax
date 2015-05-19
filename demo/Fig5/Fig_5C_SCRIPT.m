% script to rank n-grams by "informativeness" for comparisons between wild
% isolates (mutual information or F-statistic) and to see how
% informativeness relates to rank
%
% Reproduces Fig. 5C

% should repeats be included?
repeats = false;

% should approximate matches be counted?
approxMatch = false;
r2Thresh = 0.8; % threshold for counting an approximate match

% set the root directory
directory = ['/Users/abrown/Andre/wormVideos/results-12-05-10/' ...
    'wild-isolates/'];

% load the representative postures
load(['/Users/abrown/Andre/wormVideos/results-12-05-10/Laura Grundy/'...
    'gene_NA/allele_NA/N2/on_food/XX/30m_wait/'...
    'postures_90-centers_20-files_5000-framesPerFile.mat'])
postures = postures';

% get the correlation coefficient for each posture pair
rMat = vectorCor(postures, postures);


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
[uniqueNames, ~, classInts] = unique(wormNames);

% size of n-gram to use
for n = 1:3
    
    % initialisation and model set up
    seqCell = cell(numel(uniqueNames), 1); % cell for storing state sequences
    nGramCell = cell(numel(uniqueNames), 1); % cell for storing n-grams
    freqCell = cell(numel(uniqueNames), 1); % frequency distributions for each strain
    uniqueNGrams = [];
    
    % set the colour matrix
    colVec = rand(numel(uniqueNames), 3);
    colVec = colVec ./ (repmat(sum(colVec, 2), 1, 3) * 0.8);
    colVec(colVec > 1) = 1;
    
    
    % find the total set of unique n-grams present in all the wild isolate data
    % and make a cell containing the different wild isolate state sequences
    for jj = 1:numel(uniqueNames)
        disp(jj/numel(uniqueNames))
        
        % get the file names of the current strain
        currentInds = find(strcmp(wormNames, uniqueNames{jj}));
        
        % loop through files to populate seqCell and get unique n-grams
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
            seqCell{jj, ii} = expandedSeq(:, 1)';
            
            % get n-grams
            nGrams = n_gramsNumerical(expandedSeq(:, 1)', n);
            nGramCell{jj, ii} = nGrams;
            
            % get the counts
            [uniqueRows, ~] = countUniqueRows(nGrams);
            
            % append current unique rows to total uniqueNGrams
            currentNGrams = [uniqueNGrams; uniqueRows];
            
            % re-calculate unique n-grams from larger set
            [uniqueNGrams, ~] = countUniqueRows(currentNGrams);
        end
    end
    
    
    % make a matrix of counts for each of the unique n-grams in each individual
    % worm's postural sequence
    countMat = NaN(numel(wormNames), size(uniqueNGrams, 1));
    freqMat = NaN(numel(wormNames), size(uniqueNGrams, 1));
    nn = 1;
    % loop over strains
    for ii = 1:numel(uniqueNames)
        disp(ii/numel(uniqueNames))
        % loop over non-empty state sequences for current strain
        for jj = 1:sum(~cellfun(@isempty, seqCell(ii, :)))
            % loop over unique n-grams
            for kk = 1:size(uniqueNGrams, 1)
                % count the number of occurences of the current unique n-gram
                % in the current state sequence
                if approxMatch
                    matchInds = findApproxFast(nGramCell{ii, jj}, ...
                        uniqueNGrams(kk, :), r2Thresh, rMat);
                    countMat(nn, kk) = sum(matchInds);
                    freqMat(nn, kk) = ...
                        sum(matchInds) / size(nGramCell{ii, jj}, 1);
                else
                    counts = ...
                        numel( strfind(seqCell{ii, jj}, uniqueNGrams(kk, :)) );
                    countMat(nn, kk) = counts;
                    freqMat(nn, kk) = counts / numel(seqCell{ii, jj});
                end
            end
            nn = nn + 1;
        end
    end
    
    
    
    % do pair-wise comparisons of strains
    pValCell = cell(numel(uniqueNames), numel(uniqueNames));
    
    % numHits = 100; % how many sequences should be kept from each comparison
    for ii = 1:numel(uniqueNames)-1
        % get the n-grams of the ii'th strain
        nGrams1 = vertcat(nGramCell{ii, :});
        
        % get name matches of ii'th strain
        nameMatches1 = find(strcmp(wormNames, uniqueNames{ii}));
        
        for jj = ii+1:numel(uniqueNames)
            disp([ii/numel(uniqueNames), jj/numel(uniqueNames)])
            % get the n-grams of the jj'th strain
            nGrams2 = vertcat(nGramCell{jj, :});
            
            % get the unique n-grams of the current pair
            [uniqueNGramsPair, ~] = countUniqueRows([nGrams1; nGrams2]);
            
            % find the indices of the current unique n-grams in the total set
            % of unique n-grams.
            [~, matchInds] = ismember(uniqueNGramsPair, uniqueNGrams, 'rows');
            
            % get the submatrix from count mat for the current pair
            nameMatches2 = find(strcmp(wormNames, uniqueNames{jj}));
            rowInds = [nameMatches1; nameMatches2];
            countMatPair = countMat(rowInds, matchInds);
            freqMatPair = freqMat(rowInds, matchInds);
            
            % drop columns corresponding to very rare behaviours
            dropInds = sum(countMatPair) < 5;
            countMatPair(:, dropInds) = [];
            freqMatPair(:, dropInds) = [];
            matchInds(dropInds) = [];
            
            % calculate the F-stats for each behaviour
            classIntsPair = [ones(numel(nameMatches1), 1); ...
                ones(numel(nameMatches2), 1) + 1];
            
            pVals = NaN(size(freqMatPair, 2), 1);
            for kk = 1:size(freqMatPair, 2)
                % also do rank sum tests
                p = ranksum(freqMatPair(classIntsPair == 1, kk), ...
                    freqMatPair(classIntsPair == 2, kk));
                
                % get the p-value
                pVals(kk) = p;
            end
            
            % add to cell
            pValCell{ii, jj} = pVals;
        end
    end
    
    % need to save and re-load later because of required memory for larger
    % values of n if there are a large number of strains to compare
    save(['Fig_5C_WORKSPACE_' num2str(n) '-gram.mat'])
end

% clear all previous calculations to be safe
clear all

% load p-values and correct for multiple comparisons
pValsTotal = [];
for n = 3
    load(['Fig_5C_WORKSPACE_' num2str(n) '-gram.mat'], 'pValCell')
    for ii = 1:size(pValCell, 1)-1
        for jj = ii+1:size(pValCell, 2)
            pValsTotal = [pValsTotal; pValCell{ii, jj}];
        end
    end
end
[~, pCrit, bhFDR] = fdr_bh(pValsTotal, 0.05 , 'dep');



% go back through each comparison and find hits based on the globally
% controlled FDR using pCrit
for n = 3
    % load the entire workspace since it will be re-saved after updated
    % below
    load(['Fig_5C_WORKSPACE_' num2str(n) '-gram.mat'])
    
    comparisonStats = cell(numel(uniqueNames), numel(uniqueNames));
    informativeBehaviours = cell(numel(uniqueNames), numel(uniqueNames));
    for ii = 1:numel(uniqueNames)-1
        % get the n-grams of the ii'th strain
        nGrams1 = vertcat(nGramCell{ii, :});
        
        % get name matches of ii'th strain
        nameMatches1 = find(strcmp(wormNames, uniqueNames{ii}));
        
        for jj = ii+1:numel(uniqueNames)
            disp([ii/numel(uniqueNames), jj/numel(uniqueNames)])
            % get the n-grams of the jj'th strain
            nGrams2 = vertcat(nGramCell{jj, :});
            
            % get the unique n-grams of the current pair
            [uniqueNGramsPair, ~] = countUniqueRows([nGrams1; nGrams2]);
            
            % find the indices of the current unique n-grams in the total set
            % of unique n-grams.
            [~, matchInds] = ismember(uniqueNGramsPair, uniqueNGrams, 'rows');
            
            % get the submatrix from count mat for the current pair
            nameMatches2 = find(strcmp(wormNames, uniqueNames{jj}));
            rowInds = [nameMatches1; nameMatches2];
            countMatPair = countMat(rowInds, matchInds);
            freqMatPair = freqMat(rowInds, matchInds);
            
            % drop columns corresponding to very rare behaviours
            dropInds = sum(countMatPair) < 5;
            countMatPair(:, dropInds) = [];
            freqMatPair(:, dropInds) = [];
            matchInds(dropInds) = [];
            
            % calculate the F-stats for each behaviour
            classIntsPair = [ones(numel(nameMatches1), 1); ...
                ones(numel(nameMatches2), 1) + 1];
            
            
            % control false discovery rate
            hitInds = find(pValCell{ii, jj} < pCrit);
            
            % sort hitInds by bhFDR
            [~, hitSortInds] = sort(pValCell{ii, jj}(hitInds));
            hitInds = hitInds(hitSortInds);
            
            % get the frequencies for each strain
            freq1 = mean(freqMatPair(classIntsPair == 1, :));
            [~, freqSortInds1] = sort(freq1, 'descend');
            freq2 = mean(freqMatPair(classIntsPair == 2, :));
            [~, freqSortInds2] = sort(freq2, 'descend');
            
            % for the hits, record the following:
            % frequency in strain 1, frequency in strain 2, rank in
            % strain 1, and rank in strain 2
            comparisonMat = NaN(length(hitInds), 5);
            
            % add data
            comparisonMat(:, 1) = freq1(hitInds)'; % frequencies in strain 1
            comparisonMat(:, 2) = freq2(hitInds)'; % frequencies in strain 2
            [~, freqRank1] = ismember(hitInds, freqSortInds1);
            comparisonMat(:, 3) = freqRank1; % rank in strain 1
            [~, freqRank2] = ismember(hitInds, freqSortInds2);
            comparisonMat(:, 4) = freqRank2; % rank in strain 2
            comparisonMat(:, 5) = pValCell{ii, jj}(hitInds); % q-values for each comparison
            
            
            % add to cell
            comparisonStats{ii, jj} = comparisonMat;
            
            % also record informative n-grams
            informativeBehaviours{ii, jj} = ...
                uniqueNGrams(matchInds(hitInds), :);
        end
    end
    
    % save the updated workspace
    save(['Fig_5C_WORKSPACE_' num2str(n) '-gram.mat'])
end
