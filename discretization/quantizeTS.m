function [stateSequence, distVec] = ...
    quantizeTS(timeSeries, centers, warpNum, flickThresh)

% QUANTIZETS Reduce a multidimensional time series to a quantized
% representation using a set of representative posture centers.  Each time
% point in the time series is replaced by the index in postureCenters of
% the corresponding closest representative.
%
% Input
%   timeSeries    - a numFrames x numDimensions matrix
%   centers       - a matrix of cluster centers used to represent each
%                   frame of timeSeries
%   warpNum       - the maximum length of repeated values allowed in
%                   stateSequence. If warpNum were set to 4, then a
%                   sequence 'aaaabbbb' would become 'ab'.  If warpNum were
%                   set to 3 it would become 'aabb', if 2 'aabb'.
%                   warpNum == 1 corresponds to no warping.
%   flickThresh   - the threshold for considering a state change to be a
%                   noisy "flicker".  A state change A -> B -> A is
%                   considered a flicker if B lasts for flickThresh or
%                   fewer frames.  In this case, B will be set to A to
%                   remove the flicker.  flickThresh == 0 indicates no
%                   flicker removal.  N.B. Using flickThresh > 1 might not
%                   make sense in some circumstances.  For example, it
%                   probably doesn't make sense to set B to A in a sequence
%                   like XXXABBBAXXX.
%
% Output
%   stateSequence - A numStates by 2 matrix.  The first column is the 
%                   sequence of states defined by the binned values
%                   across the channels, the second column is the number of
%                   times that state appears in the discrete version of the
%                   time series before warping by warpNum.
%   distVec       - The numFrames by 1 vector with the distance between
%                   each element of timeSeries and its nearest neighbour in
%                   centers. N.B. distVec is the distance calculated BEFORE
%                   removing 'flickers'.
% 
% Copyright Medical Research Council 2013
% André Brown, andre.brown@csc.mrc.ac.uk, aexbrown@gmail.com
% 
% 
% The MIT License
% 
% Copyright (c)  Medical Research Council 2015
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.



% quantize the time series
[stateSequence, distVec] = knnsearch(timeSeries, centers, 1);

% get sign of state changes
stateDiff = sign(diff(stateSequence));

% find and remove different length flickers in turn
for ii = 1:flickThresh
    % find candidate flickers
    flicks = strfind(stateDiff', [1 zeros(1, ii - 1) -1]);
    flicks = [flicks strfind(stateDiff', [-1 zeros(1, ii - 1) 1])];
    
    % check if states at beginning and end of candidates are the same
    for jj = 1:numel(flicks)
        if stateSequence(flicks(jj)) == ...
                stateSequence(flicks(jj) + ii + 1)
            % we have a flicker, mark it for removal
            stateSequence(flicks(jj) + 1:flicks(jj) + ii) = ...
                repmat(stateSequence(flicks(jj)), ii, 1);
        end
    end
    
end

% find places where the state changes
changeInds = [0; find(diff(stateSequence) ~= 0); numel(stateSequence)] + 1;

% get the lengths of the constant segments
constLengths = diff(changeInds);

% find how many repeats must be removed from each constant segment
targetLengths = ceil(constLengths/warpNum);
cutLengths = constLengths - targetLengths;

% loop through constant state segments
dropInds = zeros(size(stateSequence));
for ii = 1:numel(changeInds) - 1
    % check that this segment must be reduced
    if cutLengths(ii) > 0
        % set dropInds to 1 where constant segments must be shortened
        dropInds(changeInds(ii):changeInds(ii) + cutLengths(ii) - 1) = 1;
    end
end

% reduce repeats in stateSequence
stateSequence(dropInds == 1) = [];

% get the number of times each state in the reduced state sequence was 
% present in the original
stateLengths = NaN(size(stateSequence));
count = 1;
for ii = 1:numel(targetLengths)
    % if the target length is 1, then the state length is simply the
    % corresponding cutLength + 1
    if targetLengths(ii) == 1
        currentLengths = cutLengths(ii) + 1;
    else
        % if the target length is greater than 1, then all lengths must be
        % warpNum, except possibly the last
        modulus = mod(constLengths(ii), warpNum);
        if targetLengths(ii) > 1 && modulus == 0
            currentLengths = repmat(warpNum, targetLengths(ii), 1);
        else
            currentLengths = ...
                [repmat(warpNum, targetLengths(ii) - 1, 1); modulus];
        end
    end
    
    % add the state lengths for the current repeated segment
    stateLengths(count:count + length(currentLengths) - 1) = currentLengths;
    
    % increment count
    count = count + length(currentLengths);
end

% combine state sequence with state lengths
stateSequence = [stateSequence stateLengths];