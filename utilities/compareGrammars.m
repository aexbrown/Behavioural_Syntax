function [matchDist, overlap] = ...
    compareGrammars(gTerminals1, gTerminals2, logScale)

% COMPAREGRAMMARS quantifies the difference between sets of gTerminals. The
% distance between two grammars is the Euclidean distance between the
% indices of matches.  Non-matches are assigned the maximum index
% difference + 1 (as if the non-matched element was simply added to the end
% of the compared grammar). To weight the most compressive grammar elements
% more, the distance can also be calculated using the difference of the
% log_10 of indices.  The overlap is simply the fraction of grammar
% elements that match between each grammar.  N.B. because grammars may have
% different numbers of elements, matchDist is not symmetric, although
% overlap is.
%
% Input
%   gTerminals1,2 - a numRules by 1 cell of row vectors of terminal symbols
%                   in a grammar.
%   logScale      - logical.  If true, the index differences will be log_10
%                   transformed in the distance calculation.  This
%                   effectively increases the weight of the most
%                   compressive (low index) elements of the grammars.
%
% Output
%   matchDist     - The Euclidean distance between match indices.
%   overlap       - a matrix of overlap counts between the two sets of
%                   gTerminals.
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
indexDiffs = NaN(numel(gTerminals1), 1);
overlap = 0;

% get the lengths of the elements in gTerminals2
lengths2 = cellfun('length', gTerminals2);

% loop through the elements of the first set of grammar terminals
for ii = 1:numel(gTerminals1)
    % get the current sequence from gTerminals1
    currentSeq = gTerminals1{ii};
    
    % find elements in gTerminals2 with the same length as currentSeq
    lenMatchInds = find(numel(currentSeq) == lengths2);
    
    if any(lenMatchInds)
        % look for a match in gTerminals2
        matchInd = lenMatchInds( cellfun(@(x) all(x == currentSeq), ...
            gTerminals2(lenMatchInds)) );
        
        
        % check for multiple matches
        if numel(matchInd) > 1
            error(['multiple matches found.  '...
                'gTerminals should have only unique elements.'])
        end
        
        if ~isempty(matchInd)
            % a match was found.  Increment overlap count.
            overlap = overlap + 1;
            
            % Get the index difference between currentSeq and its match
            if logScale
                indexDiffs(ii) = log10(matchInd) - log10(ii);
            else
                indexDiffs(ii) = matchInd - ii;
            end
            
        else % no matches found
            if logScale
                indexDiffs(ii) = log10(numel(gTerminals2) + 1) - log10(ii);
            else
                indexDiffs(ii) = numel(gTerminals2) + 1 - ii;
            end
        end
        
    else % no length matches
        if logScale
            indexDiffs(ii) = log10(numel(gTerminals2) + 1) - log10(ii);
        else
            indexDiffs(ii) = numel(gTerminals2) + 1 - ii;
        end
    end
end

% calculate the distance from the index differences
matchDist = sqrt( sum(indexDiffs.^2) );
