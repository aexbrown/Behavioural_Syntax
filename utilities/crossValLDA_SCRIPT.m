% compare LDA cross-validation classification accuracy for the different
% grammars for a range of features each

% set the directory and get the list of files
directory = ['/Users/abrown/Andre/wormVideos/results-12-05-10/' ...
    'analysisResults/'];
[fileList, ~] = dirSearch(directory, 'featMat_gTerminals*.mat');

% initialise
featNumVec = [1 5 10 50 100 150 200 250 300 350 400 450 500 ...
              550 600 650 700 750 800 850 900 950 1000];
errorMat = NaN(numel(fileList), numel(featNumVec));

% load the worm names
wormNames = struct2cell(load([directory, 'wormNames.mat']));
wormNames = wormNames{:};

% loop through feature matrices
for ii = 1:numel(fileList)
    disp(ii/numel(fileList))
    % load featMat and featutreIndices
    featMat = cell2mat(struct2cell(load(fileList{ii})));
    featureIndices = cell2mat(struct2cell(load(...
        strrep(fileList{ii}, 'featMat', 'featInds'))));
    
    % loop through feature numbers
    for jj = 1:numel(featNumVec)
        % set number of features to use
        featNum = featNumVec(jj);
        
        % perform 10 fold cross validation.
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        % last class has only one member.  Leave it off for this run.
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        lda = ClassificationDiscriminant.fit( ...
            featMat(1:end-1, featureIndices(1:featNum)), wormNames(1:end-1));
        cvmodel = crossval(lda, 'kfold', 10);
        cverror = kfoldLoss(cvmodel);
        
        % add error to output
        errorMat(ii, jj) = cverror;
    end  
end