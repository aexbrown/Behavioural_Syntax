function [matchN, meanDist, stdDist, meanT, stdT] = ...
    patternMatch(stateSeq, gExp, distTable, thresh)

% PATTERNMATCH finds the number of non-overlapping instances of each
% pattern represented by sequences of terminal symbols represented in the
% expanded grammar gExp.  Approximate matches (those with distances less
% than or equal to thresh) are allowed.  The mean of standard deviation of
% these approximate matches is also output.
%
% Inputs
%   stateSeq  - A 2 by numStates matrix.  The first column is the 
%               sequence of states defined by the binned values
%               across the channels, the second column is the number of
%               times that state appears in the discrete version of the
%               time series before warping by warpNum.
%   gExp      - a number of rules by 1 cell array.  It is either the 
%               expanded version of the input grammar obtained by applying 
%               the rules encoded in the grammar or a version containing 
%               only the terminal rules.  If gExp{ii} is a vector it must
%               contain only terminal symbols.  If gExp{ii} is a cell, then
%               gExp{ii}{end} contains only terminal symbols.
%   distTable - a lookup table of distances.  distTable(i, j) is the
%               distance between state i and state j.
%   thresh    - the distance threshold for calling matches.  Distances
%               less than or equal to thresh are considered matches.
%               N.B. this is the threshold PER MATCH, so for a 5-gram,
%               the threshold will be 5*thresh
%
% Output
%   matchN   - a vector of counts.  matchN(i) contains the count of the ith
%              element in grammar.
%   meanDist - a vector indicating the mean distance between each grammar
%              rule and its approximate matches in dataVec
%   stdDist  - a vector indicating the standard deviation of the distances
%              between each grammar rule and its matches in dataVec
%   meanT    - the average time spent in the matched sequence (each time is
%              simply the sum of number of frames spent in each of the 
%              states of the matched pattern)
%   stdT     - the standard deviation of the distribution of times spent in
%              the matched sequences
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



% separate the sequence of states and the times
dataVec = stateSeq(1, :);
timeVec = stateSeq(2, :);

% initialisation
matchN = zeros(size(gExp, 1), 1);
meanDist = zeros(size(gExp, 1), 1);
stdDist = zeros(size(gExp, 1), 1);
meanT = zeros(size(gExp, 1), 1);
stdT = zeros(size(gExp, 1), 1);

% check each rule in the grammar and take most compressive
for ii = 1:size(gExp, 1)
    % get the terminal symbols in the current rule
    if iscell(gExp{ii})
        % in this case, gExp has non-terminal rules, so take the last one
        currentSeq = gExp{ii}{end};
    else
        % in this case the non-terminals have already been stripped
        currentSeq = gExp{ii};
    end
    
    nn = length(currentSeq);
    
    % get all of the n-grams
    nGrams = n_gramsNumerical(dataVec, nn);
    
    % get all the distances
    distVec = distSeq(currentSeq, nGrams, distTable);
    
    % get the matching indices
    matchInds = find(distVec <= nn * thresh);
    
    % remove the overlaps
    while any(diff(matchInds) < nn)
        % find the first overlap
        overlapInd1 = find(diff(matchInds) < nn, 1, 'first') + 1;
        matchInds(overlapInd1) = [];
    end
    
    % add match number and distances stats to output
    matchN(ii) = numel(matchInds);
    meanDist(ii) = mean(distVec(matchInds));
    stdDist(ii) = std(distVec(matchInds));
    
    % get the match times
    matchTimes = NaN(numel(matchInds), 1);
    for jj = 1:numel(matchInds)
        matchTimes(jj) = sum(timeVec(matchInds(jj):matchInds(jj)+nn-1));
    end
    meanT(ii) = mean(matchTimes);
    stdT(ii) = std(matchTimes);
end