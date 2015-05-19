function expandedSeq = expandSequence(stateSequence, dsRate)

% EXPANDSEQUENCE takes a stateSequence matrix with a sequence of states in
% the first column and a duration in the second column and returns expands
% it so that each duration is 1.  For example, the state sequence
%           [1, 2; 7, 3] becomes [1, 1; 1, 1; 7, 1; 7, 1; 7, 1]
% The function also takes into account downsampling (if any) of the
% original signal.  If the downsampling rate in the original were 2 (only
% every 2nd point included in the state sequence), then the state sequence
% [1, 2; 7, 3] becomes [1, 1, 1, 1, 7, 7, 7, 7, 7, 7]' (first column only)
% 
% Input
%   stateSequence - an nx2 matrix.  The first column is the sequence of
%                   states and the second column is the number of times
%                   they occured in the original sequence before warping.
%   dsRate        - the downsampling rate used to make the original state
%                   sequence.  A dsRate of 2 means only every second point
%                   was included.
% 
% Outpu
%   expandedSeq   - an mx2 matrix.  The unwarped version of stateSequence
%                   where every every duration is 1.
% 
% 
% André Brown, andre.brown@csc.mrc.ac.uk, aexbrown@gmail.com
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



% initialise
expandedSeq = NaN(sum(stateSequence(:, 2)), 1);

% loop through states
index = 1;
for ii = 1:size(stateSequence, 1)
    % expand the current state
    expandedState = ones(stateSequence(ii, 2), 1) * stateSequence(ii, 1);
    
    % add to expandedSeq and increment index
    expandedSeq(index:index + stateSequence(ii, 2) - 1) = expandedState;
    index = index + stateSequence(ii, 2);
end

% undo downsampling unless dsRate is 1
if dsRate ~= 1
    pointSpacing = ...
        (size(expandedSeq, 1) - 1) / (dsRate*size(expandedSeq, 1));
    expandedSeq = ...
        interp1((1:size(expandedSeq, 1))', expandedSeq(:,1), ...
        (1:pointSpacing:size(expandedSeq, 1)-pointSpacing)', 'nearest');
end

% add durations (now all ones)
expandedSeq = [expandedSeq, ones(size(expandedSeq, 1), 1)];

