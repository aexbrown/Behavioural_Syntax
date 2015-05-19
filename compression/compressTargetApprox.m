function [matchN, meanDist, stdDist, compVec] = ...
    compressTargetApprox(dataVec, grammar, distTable, thresh)

% COMPRESSTARGETAPPROX uses a grammar derived from another sequence to 
% compress a target data vector that was not used to derive the grammar. It
% calculates the compressed version of the target vector as well as the
% number of occurances of each of the elements of the grammar.
%
% See also: compressSequenceApprox.m
%
% Inputs
%   dataVec  - a row vector of numbers to be compressed
%   grammar  - a number of rules by 2 cell array. The first column has the
%              left hand side of each replacement rule while the second
%              column has the right hand side (so the first column lists
%              all non-terminals in the grammar).
%   distTable   - a lookup table of distances.  distTable(i, j) is the
%                 distance between state i and state j.
%   thresh      - the distance threshold for calling matches.  Distances
%                 less than or equal to thresh are considered matches.
%                 N.B. this is the threshold PER MATCH, so for a 5-gram,
%                 the threshold will be 5*thresh
%
% Output
%   matchN   - a vector of counts.  matchN(i) contains the count of the ith
%              element in grammar.
%   meanDist - a vector indicating the mean distance between each grammar
%              rule and its approximate matches in compVec
%   stdDist  - a vector indicating the standard deviation of the distances
%              between each grammar rule and its matches in compVec
%   compVec  - the vector that has been compressed using grammar.  dataVec
%              can be recovered by applying the grammar rules in reverse.
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



% initialisation
matchN = zeros(size(grammar, 1), 1);
meanDist = zeros(size(grammar, 1), 1);
stdDist = zeros(size(grammar, 1), 1);
compVec = dataVec;

% add enough new rows and columns to distTable to accommodate the
% the largest index in grammar.  Only exact matches are allowed for
% non-terminals so set distance above threshold for all new matches
% except self-matches which of course have zero distance.
distTable = [distTable ...
    ones(size(distTable, 1), ...
    grammar{end, 1} - size(distTable, 1)...
    ) + thresh];
distTable = [distTable; ...
    ones(grammar{end, 1} - size(distTable, 1), ...
    size(distTable, 2) ...
    ) + thresh];
for kk = grammar{end, 1}:size(distTable, 1)
    distTable(kk, kk) = 0;
end


% make all replacements until none remaining is compressive
isCompressive = true;
while isCompressive
    % reset best savings
    bestSavings = 0;
    bestCount = 0;
    sequence = [];
    bestInd = -1;
    
    % check each rule in the grammar and take most compressive
    for ii = 1:size(grammar, 1)
        % get the current sequence to check
        currentSeq = grammar{ii, 2};
        nn = length(currentSeq);
        
        % get all of the n-grams
        nGrams = n_gramsNumerical(compVec, nn);
        
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
        
        % calculate number of matches without overlap
        count = numel(matchInds);
        
        % calculate the savings of the current element
        savings = (length(currentSeq) - 1) * (count - 1) - 2;
        if savings > bestSavings
            bestSavings = savings;
            bestCount = count;
            bestMean = mean(distVec(matchInds));
            bestStd = std(distVec(matchInds));
            sequence = currentSeq;
            inds = matchInds;
            bestInd = ii; % index of best sequence in grammar
        end
    end
    
    % if no sequence was found with positive savings, compression is
    % finished
    if bestSavings == 0
        isCompressive = false;
    else
        % include the number of matches of the most compressive grammar
        % element in matchN
        matchN(bestInd) = bestCount;
        meanDist(bestInd) = bestMean;
        stdDist(bestInd) = bestStd;
        
        % make the replacements in compVec
        for jj = 1:numel(inds)
            compVec(inds(jj):inds(jj) + numel(sequence) - 1) = ...
                [grammar{bestInd, 1} NaN(1, numel(sequence) - 1)];
        end
        compVec(isnan(compVec)) = [];
    end
end