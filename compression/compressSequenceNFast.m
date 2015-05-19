function [grammar, compVec, totSavings] = ...
    compressSequenceNFast(dataVec, newStart, nMax)

% COMPRESSSEQUENCE Recursively finds the most compressive subsequence in
% dataVec and creates and replaces it with a new number.  This replacement
% creates a new rule in the grammar.  Replacements are made until there are
% none left that lead to further compression.  See the following paper
% for more details: Nevill-Manning and Witten (2000) On-Line and Off-Line
% Heuristics for Inferring Hierarchies of Repetitions in Sequences.
% Proceedings of the IEEE 88:1745.
%
% Input
%   dataVec  - a row vector of numbers to be compressed
%   newStart - this is the number that will be used to label the first new
%              rule in the grammar.  It must be greater than the maximum
%              value in dataVec.  If empty, then max(dataVec) + 1 is used.
%   nMax        - the maximum length n-gram to check for compression
%
% Output
%   grammar  - a number of rules by 2 cell array. The first column has the
%              left hand side of each replacement rule while the second
%              column has the right hand side (so the first column lists
%              all non-terminals in the grammar).
%   compVec  - the vector that has been compressed using grammar.  dataVec
%              can be recovered by applying the grammar rules in reverse.
%   totSavings - the total space saving achieved during the compression,
%                taking into account the size of the created grammar rules
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



% check dataVec
if size(dataVec, 1) ~= 1
    error('dataVec must be a row vector.')
end

% define newStart if left empty
if isempty(newStart)
    newStart = max(dataVec) + 1;
end

% check that newStart is large enough
if newStart <= max(dataVec)
    error('newStart must be greater than max(dataVec).')
end

% initialise grammar
grammar = cell(1, 2);

% initialise compVec and make a suffix array
compVec = dataVec;
totSavings = 0;

% compress segments until none are found that lead to compression
sequence = NaN;
newInd = newStart;
ii = 1;
while ~isempty(sequence)
    % find the most compressive sequence in dataVec
    [sequence, locations, savings] = compressiveNFast(compVec, nMax);
    
    % update the total savings (i.e. compression)
    totSavings = totSavings + savings;
    
    % add the rule to grammar
    grammar{ii, 1} = newInd;
    grammar{ii, 2} = sequence;
    ii = ii + 1;
    % disp(ii);
    
    % make the replacements.  Note: strrep does not work here.  For example
    % if sequence is [44 68 44] and compVec has a subsequence that is 
    % [44 68 44 68 44 68 44 448], strrep will give [68 480 480 480 448]
    % which is wrong.
    for jj = 1:numel(locations)
        compVec(locations(jj):locations(jj) + numel(sequence) - 1) = ...
            [newInd NaN(1, numel(sequence) - 1)];
    end
    compVec(isnan(compVec)) = [];
    
    newInd = newInd + 1;
end

% remove the last (empty) entry of the grammar
grammar = grammar(1:end-1, :);
