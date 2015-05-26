function [sequence, locations, bestSavings] = ...
    compressiveNFast(dataVec, nMax)

% COMPRESSIVEN Finds the sequence in dataVec that results in the largest
% compression when it is replaced by a new rule, but only searches n-grams
% up to nMax to save time.  See the following paper for more details on the
% non-approximate version:
% Nevill-Manning and Witten (2000) On-Line and Off-Line
% Heuristics for Inferring Hierarchies of Repetitions in Sequences.
% Proceedings of the IEEE 88:1745.
%
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% N.B. In the case of ties, where two subsequences would give equal
% compression, Nevill-Manning takes the one that occurs first in the
% sequence, whereas here we take the one that occurs first in the sorted
% unique n-grams.  This can change the grammar that is produced by
% compressSequenceNFast and the total compression achieved.
%
% There may also be similar differences between using compressiveN and
% compressiveNFast because of the way unique() sorts strings vs. numerical
% values.
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%
% SEE ALSO: compressive.m
%
% Input
%   dataVec     - a row vector of numbers to be compressed
%   nMax        - the maximum length n-gram to check for compression
%
% Output
%   sequence    - the sequence in dataVec that results in the best
%                 compression when replaced with a new rule (best
%                 compression balances a high frequency against a long
%                 sequence)
%   locations   - the locations where sequence occurs
%   bestSavings - the compression achieved by making the replacements (this
%                 takes into account the size of the rule that is created)
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


% check that dataVec is not too short.  If it is, only search for n-grams
% up to length of dataVec - 1
if length(dataVec) <= nMax
    disp(['Warning: nMax was reduced to ' num2str(length(dataVec) - 1)])
    nMax = length(dataVec) - 1;
end

% initialise
bestSavings = 0;
sequence = [];

% convert dataVec to a string with double spaces between entries
dataString = ['  ' sprintf('%1d  ', dataVec)];

% loop through n-grams
for nn = 2:nMax
    % get all of the n-grams
    nGrams = n_gramsNumerical(dataVec, nn);
    
    % get the count of the unique n-grams.  This count includes overlaps,
    % but runs much faster than the regexp search for long strings.  We
    % check for overlaps below
    [uniqueNGrams, countOverlaps] = countUniqueRows(nGrams);
    
    % check if there is only 1 unique n-gram (a repeated sequence
    if size(uniqueNGrams, 1) == 1
        nGramString = [' ' sprintf('%1d  ', uniqueNGrams)];
        nGramString = nGramString(1:end-1);
        bestCount = numel( regexp(dataString, nGramString) );
        bestNGram = uniqueNGrams;
    else
        
        % find the most frequently occuring n-grams
        [~, sortInds] = sort(countOverlaps, 'descend');
        
        % counts without overlaps are less than or equal to counts with
        % overlaps for any given sequence.  So check the top hits with overlaps
        % and compare the counts without overlaps.  If the count without
        % overlaps is >= to the next best count _with_ overlaps, it is the
        % best.  If not, check the next-best hit with overlaps for its count
        % without overlaps.
        
        previousBest = 0;
        previousBestInd = 1;
        for ii = 1:numel(sortInds)
            % if the current count with overlaps is no better than the previous
            % best count without overlaps, then stop and take the previous best
            if previousBest >= countOverlaps(sortInds(ii))
                bestNGram = uniqueNGrams(sortInds(previousBestInd), :);
                bestCount = previousBest;
                break
            end
            
            % get the non-overlapping counts of the current n-gram. Pad
            % uniqueNGrams with spaces to avoid matches like: '76 7' matching
            % '76 78' in dataString.  Each element of dataVec is separated
            % by two spaces, so separate each element of the current n-gram
            % with two spaces but only pad with one on either side.  This
            % allows overlaps to be correctly counted.  E.g. ' 2  2 ' will
            % be found twice in '12  2  2  2  2  22  5' as it should be.
            nGramString = [' ' sprintf('%1d  ', uniqueNGrams(sortInds(ii), :))];
            nGramString = nGramString(1:end-1);
            countNoOverlaps = numel( regexp(dataString, nGramString) );
            
            % check if you have gotten to the last possible n-gram
            if ii == numel(sortInds)
                % this means none of the previous n-grams was best.  This is
                % unlikely on real data, but can happen if there is a tie in
                % the counts across all n-grams.
                bestNGram = uniqueNGrams(sortInds(previousBestInd), :);
                bestCount = previousBest;
                break
            end
            
            % check if the counts without overlaps of the current n-gram is the
            % highest of any n-gram.
            if countNoOverlaps >= countOverlaps(sortInds(ii + 1))
                % this is the most frequent n-gram
                bestNGram = uniqueNGrams(sortInds(ii), :);
                bestCount = countNoOverlaps;
                break
                % check if it is better than the previous best
            elseif countNoOverlaps > previousBest
                previousBest = countNoOverlaps;
                previousBestInd = ii;
            end
        end
        
        % for testing, check all unique N-grams and take best over all.
        % Make sure it's the same as that found above
        %     nGramsTest = n_grams(dataString, nn);
        %     [uniqueNGramsTest, ~] = countUniqueStr(nGramsTest);
        %     counts = NaN(numel(uniqueNGramsTest), 1);
        %     for ii = 1:numel(uniqueNGramsTest)
        %         counts(ii) = numel( regexp(dataString, ...
        %             [' ' uniqueNGramsTest{ii} ' ']) );
        %     end
        %     if max(counts) ~= bestCount
        %         error('brute force does not agree with early abandon solution')
        %     end
        
    end
    % calculate the savings of using the current best n-gram for
    % compression
    savings = (nn - 1) * (bestCount - 1) - 2;
    if savings > bestSavings
        bestSavings = savings;
        sequence = bestNGram;
    end
    
end

% get the locations with overlaps
locations = strfind(dataVec, sequence);

% remove the overlaps
while any(diff(locations) < numel(sequence))
    % find the first overlap
    overlapInd1 = find(diff(locations) < numel(sequence), 1, 'first') + 1;
    locations(overlapInd1) = [];
end