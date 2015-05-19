function [matchN, compVec] = compressTarget(dataVec, grammar)

% COMPRESSTARGET uses a grammar derived from another sequence to compress a
% target data vector that was not used to derive the grammar.  It
% calculates the compressed version of the target vector as well as the
% number of occurances of each of the elements of the grammar.
%
% See also: compressSequenceNFast.m and expandGrammar.m
%
% Inputs
%   dataVec  - a row vector of numbers to be compressed
%   grammar  - a number of rules by 2 cell array. The first column has the
%              left hand side of each replacement rule while the second
%              column has the right hand side (so the first column lists
%              all non-terminals in the grammar).
%
% Output
%   matchN   - a vector of counts.  matchN(i) contains the count of the ith
%              element in grammar.
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
compVec = dataVec;

% make all replacements until none remaining is compressive
isCompressive = true;
while isCompressive
    % reset best savings
    bestSavings = 0;
    bestCount = 0;
    sequence = [];
    bestInd = -1;
    
    % make a string version of compVec for the search in the loop
    dataString = [' ' sprintf('%1d ', compVec)];
    
    % check each rule in the grammar and take most compressive
    for ii = 1:size(grammar, 1)
        % get the current sequence to check
        currentSequence = grammar{ii, 2};
        
        % first do fast check to see if the sequence is even present before
        % doing slower no-overlap  count
        if ~isempty( strfind(compVec, currentSequence) )
            
            % get the non-overlapping counts of the current n-gram. Pad
            % uniqueNGrams with spaces to avoid matches like: '76 7'
            % matching '76 78' in dataString.
            currentString = [' ' sprintf('%1d ', currentSequence)];
            count = numel( regexp(dataString, currentString) );
            
            % calculate the savings of the current element
            savings = (length(currentSequence) - 1) * (count - 1) - 2;
            if savings > bestSavings
                bestSavings = savings;
                bestCount = count;
                sequence = currentSequence;
                bestInd = ii; % index of best sequence in grammar
            end
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
        
        % get the locations with overlaps
        inds = strfind(compVec, sequence);
        
        % remove the overlaps
        while any(diff(inds) < numel(sequence))
            % find the first overlap
            overlapInd1 = ...
                find(diff(inds) < numel(sequence), 1, 'first') + 1;
            inds(overlapInd1) = [];
        end
        
        % make the replacements in compVec
        for jj = 1:numel(inds)
            compVec(inds(jj):inds(jj) + numel(sequence) - 1) = ...
                [grammar{bestInd, 1} NaN(1, numel(sequence) - 1)];
        end
        compVec(isnan(compVec)) = [];
    end
end