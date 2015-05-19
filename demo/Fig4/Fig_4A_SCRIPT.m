% make a plot of worm trajectories in different conditions.  Paths are
% oriented along their first principal component and started from (0, 0).
% 
% Reproduces Fig. 4A


% set the root directory where the projected amplitudes are
directory = '/Users/abrown/Andre/wormVideos/robynTracking/no-food/';
% directory = '/Users/abrown/Andre/wormVideos/robynTracking/benzaldehyde/';
% directory = ['/Users/abrown/Andre/wormVideos/results-12-05-10/Laura Grundy/' ...
%     'gene_NA/allele_NA/N2/on_food/XX/30m_wait/'];

% get file names
[fileList, ~] = dirSearch(directory, 'features.mat');

for ii = 1:numel(fileList)
    load(fileList{ii})
    
    % load x and interpolate over NaNs
    x = worm.path.coordinates.x;
    x(isnan(x)) = interp1(find(~isnan(x)),...
        x(~isnan(x)), find(isnan(x)),'linear');
    
    % drop any remaining NaNs
    x(isnan(x)) = [];
    
    % load y and interpolate over NaNs
    y = worm.path.coordinates.y;
    y(isnan(y)) = interp1(find(~isnan(y)),...
        y(~isnan(y)), find(isnan(y)), 'linear');
    
    % drop any remaining NaNs
    y(isnan(y)) = [];
    
    if any(~isnan(x))
        
        % get the principal components of the path
        coeffs = pca([x; y]');
        
        % make a rotation matrix from the angle of the first PC
        theta = -atan(coeffs(1, 1)/coeffs(2, 1));
        rotationMat = [cos(theta), -sin(theta); ...
            sin(theta), cos(theta)];
        
        % rotate trajectory
        rotCoords = [x; y]'*rotationMat;
        
        % there is a sign ambiguity for the PCs.  If the trajectory ends up
        % pointing down, rotate by pi
        if rotCoords(1, 2) > rotCoords(end, 2)
            rotCoords = rotCoords * [-1, 0; 0, -1];
        end
        
        % plot trajectory starting with first point at zero
        line(rotCoords(:, 1) - rotCoords(1, 1), ...
            rotCoords(:, 2) - rotCoords(1, 2), 'LineWidth', 2)
        %     xlim([-20000 20000])
        ylim([-30000 40000])
        axis equal
    end
end

% add a dot at zero to indicate start
line(0, 0, 'MarkerSize', 20, 'Marker', '.', 'Color', 'r')
