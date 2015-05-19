% cluster a matrix of skeleton angles into a fixed number of templates
% using k-means

% load representative postures
load(['/Users/abrown/Andre/wormVideos/results-12-05-10/Laura Grundy/' ...
    'gene_NA/allele_NA/N2/on_food/XX/30m_wait/' ...
    'postures_90-centers_20-files_5000-framesPerFile.mat'])
postures = postures';

% how much should the data be downsampled?
ds = 5;

% should skeletons be rotated to reflect original data or plotted with
% constant angle?
rotateSkeleton = 1;


%--------------------------------------------------------------------------
% load the data
%--------------------------------------------------------------------------

% set the root directory
directory = ['/Users/abrown/Andre/wormVideos/results-12-05-10/'...
    'Laura Grundy/gene_NA/allele_NA/N2/on_food/XX/30m_wait/'];

% get the list of files (includes data from mutant classes)
[fileList, ~] = dirSearch(directory, 'angleArray.mat');

% choose file to plot from list
fileNum = 10;

% load data
angleArray = cell2mat(struct2cell(load(fileList{fileNum}, 'angleArray')));
meanArray = cell2mat(struct2cell(load(...
    strrep(fileList{fileNum}, 'angleArray', 'meanAngles'), 'meanAngles')));
angleArray = angleArray(1:ds:end, :);
meanArray = meanArray(1:ds:end);

% if the worm is in the "left" folder, invert all the angles.
if ~isempty(strfind(fileList{fileNum}, '/L/'))
    angleArray = angleArray * -1;
end


% instead of dropping NaN values, interpolate over NaNs
% initialise anglesNoNaN
anglesNoNaN = angleArray;

% interpolate over NaN values
for jj = 1:size(angleArray, 2)
    pAmp = angleArray(:, jj);
    pAmp(isnan(pAmp)) = interp1(find(~isnan(pAmp)),...
        pAmp(~isnan(pAmp)), find(isnan(pAmp)),'linear', 'extrap');
    anglesNoNaN(:, jj) = pAmp;
end

% fix direction jumps.  Use nearest neighbour interpolation to find jumps
meansTemp = meanArray;
meansTemp(isnan(meansTemp)) = interp1(find(~isnan(meansTemp)),...
    meansTemp(~isnan(meansTemp)), find(isnan(meansTemp)), ...
    'nearest', 'extrap');
jumpInds = find(abs(diff(meansTemp)) > 1.5*pi);
jumpInds = jumpInds + 1;

% correct jumps
meansNoNaN = meanArray;
for ii = 1:length(jumpInds)
    meansNoNaN(jumpInds(ii):end) = ...
        meansNoNaN(jumpInds(ii):end) - ...
        sign(meansTemp(jumpInds(ii)) - meansTemp(jumpInds(ii)-1)) * 2*pi;
end

% now interpolate linearly to get smooth change during e.g. omegas
meansNoNaN(isnan(meansNoNaN)) = interp1(find(~isnan(meansNoNaN)),...
    meansNoNaN(~isnan(meansNoNaN)), find(isnan(meansNoNaN)), ...
    'linear', 'extrap');



%--------------------------------------------------------------------------
% export image series showing the fit quality for various values of k
%--------------------------------------------------------------------------



% % make a movie some skeletons over time with their nearest neighbours
% vidObj = VideoWriter(['postures_' num2str(k(kk)) '.avi']);
% vidObj.FrameRate = 20;
% vidObj.Quality = 100;
%
% open(vidObj);

scaleFactor = 1.5;
figure('Position',[100 100 740 581]*scaleFactor)
% figure
for ii = 100:400
    % get the matching center
    [centerInd, ~] = knnsearch(anglesNoNaN(ii, :), postures, 1);
    
    % calculate skeleton reconstructions
    if rotateSkeleton
        [x, y] = angle2skel(anglesNoNaN(ii, :)', -meansNoNaN(ii), 1.7);
        [xFit, yFit] = ...
            angle2skel(postures(centerInd, :)', -meansNoNaN(ii), 1.7);
    else
        [x, y] = angle2skel(anglesNoNaN(ii, :)', 0, 1.7);
        [xFit, yFit] = ...
            angle2skel(postures(centerInd, :)', 0, 1.7);
    end
    
    % also plot the corresponding worm skeleton on the left
    plot((xFit - mean(xFit)), yFit - mean(yFit), '.', ...
        'Color', [0.2 0.4 1], 'MarkerSize', 40 * scaleFactor)
    hold on
    % mark the head
    plot((xFit(1) - mean(xFit)) , ...
        yFit(1) - mean(yFit),...
        'o', 'Color', [0.2 0.4 1], 'MarkerSize', 30 * scaleFactor)
    plot((xFit(1) - mean(xFit)) , ...
        yFit(1) - mean(yFit),...
        'o', 'Color', [0.2 0.4 1], 'MarkerSize', 21 * scaleFactor)
    
    plot((x - mean(x)), y - mean(y), '.', ...
        'Color', [1 0.1 0.1], 'MarkerSize', 40 * scaleFactor)
    % mark the head
    plot((x(1) - mean(x)), y(1) - mean(y), 'o',...
        'Color', [1 0.1 0.1], 'MarkerSize', 30 * scaleFactor)
    plot((x(1) - mean(x)), y(1) - mean(y), 'o',...
        'Color', [1 0.1 0.1], 'MarkerSize', 21 * scaleFactor)
    
    xlim([-1, 1])
    ylim([-1, 1])
    axis equal
    
    set(gcf, 'Color', 'w')
    axis off
    
    %     pause(0.01)
    %     F = capturescreen(gcf);
    %     writeVideo(vidObj,F);
    
    saveas(1, ...
        ['postures-' num2str(90) '_ds-' num2str(ds) ...
        '_frame-' num2str(ii) '.png'])
    hold off
end
% close(vidObj);
close all