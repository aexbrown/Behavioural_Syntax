function [matchInds1, matchInds2] = partialMatchRow(mat1, mat2)

% PARTIALMATCHROW checks whether any of the rows of mat2 matches any of the
% sub-rows of mat1 (or entire rows if mat1 and mat2 have the same number of
% columns)
% 
% mat1: 10    81    65  
%       40    65    76
%       42    24    50
%       30    64    19
%       64    19    46
%
% mat2: 24    50
%       30    64
%       -1    -1
%
% returns matchInds: [0; 0; 1; 1; 0]
%         unMatchedInds: [1; 1; 0]
%
% Input
%   mat1, 2       - matrices to check for row matches
% 
% Output
%   matchInds     - logical vector indicating any matching rows in mat1
%   unMatchedInds - logical vector indicating any unmatched rows in mat2
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



% get the number of columns
numCols1 = size(mat1, 2);
numCols2 = size(mat2, 2);

% check inputs
if numCols1 < numCols2
    error('mat1 must have more columns than mat2')
end

% initialise
matchInds1 = zeros(size(mat1, 1), 1);
matchInds2 = zeros(size(mat2, 1), 1);

% loop through sub-rows of mat1
for ii = 1:numCols1 - numCols2 + 1
    % get current partial matches
    [currentMatchInds, currentUnMatchedInds] ...
        = ismember(mat1(:, ii:ii+numCols2-1), mat2, 'rows');
    
    % update over all matchInds
    matchInds1 = matchInds1 | currentMatchInds;
    
    % non-zero elements of the 2nd output of ismember contain the locations
    % in mat2 that were matched
    mat2MatchInds = ...
        unique(currentUnMatchedInds(currentUnMatchedInds ~= 0));
    matchInds2(mat2MatchInds) = 1;
end
    