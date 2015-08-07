function [skelX, skelY] = angle2skel(angleArray, meanAngle, arclength)

% ANGLE2SKEL Take in an angle array and integrate over the angles to get
% back a skeleton.  NB: This reconstruction assumes each segment was
% equally spaced so that each reconstructed skeleton segment has length
% arclength/(numAngles + 1)
% 
% Input
%   angleArray - a numSkelPoints - 1 by numFrames array of skeleton
%                tangent angles that have been rotated to have a mean
%                angle of zero.
%   meanAngle  - a 1 by numFrames array of angles.  Each angle is the mean
%                angle of the skeleton used to make the corresponding row
%                of angles in angleArray.  Can be left as zeros if no
%                rotation is desired.
%   arclength  - the total arclength of the skeleton to be reconstructed.
%                Can be set to 1 for a normalised skeleton.
% 
% Output
%   skelX      - a numAngles + 1 by numFrames array of skeleton
%                x-coordinates
%   skelY      - a numAngles + 1 by numFrames array of skeleton
%                y-coordinates

% get dimensions
numAngles = size(angleArray, 1);
numFrames = size(angleArray, 2);


% initialisation
skelX = NaN(numAngles + 1, numFrames);
skelY = NaN(numAngles + 1, numFrames);

for ii = 1:numFrames
    % add up x-contributions of angleArray, rotated by meanAngle
    skelX(:, ii) = [0; cumsum(cos( angleArray(:, ii) + meanAngle(ii) ) * ...
        arclength/numAngles) ];
    
    % add up y-contributions of angleArray, rotated by meanAngle
    skelY(:, ii) = [0; cumsum(sin( angleArray(:, ii) + meanAngle(ii) ) * ...
        arclength/numAngles) ];
end

% % calculate actual acrlength
% segLength = sum(sqrt(diff(skelX).^2 + diff(skelY).^2));
% disp(segLength)
% disp(arclength)

