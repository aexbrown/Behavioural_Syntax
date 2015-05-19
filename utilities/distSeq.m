function dist = distSeq(seq1, nGrams, distTable)

% DISTSEQ finds the distsance between two discrete distances using a lookup
% table of distances between each element in the sequence
% 
% Input
%   seq1      - a row vector of integers.  Each state in the sequence is an
%               index for looking up distances in distTable
%   nGrams    - a number of sequences x n matrix of n-grams.
%               size(nGrams, 2) must equal size(seq1, 2)
%   distTable - a lookup table of distances.  distTable(i, j) is the
%               distance between state i and state j.
% 
% Output
%   dist      - a vector of distances between seq1 and each row of nGrams
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



% check inputs
if size(seq1, 2) ~= size(nGrams, 2)
    error('input sequence must have same length as n-grams in nGram matrix')
end

% get the linear indices in distTable defined by the sequences (bsxfun
% version is significantly faster than the repmat/sub2ind version commented
% out below)
seqInds = bsxfun(@plus, size(distTable, 1) * (seq1 - 1), nGrams);

% seqInds = ...
%     sub2ind(size(distTable), repmat(seq1, size(nGrams, 1), 1), nGrams);

% sum distance between seq1 and nGrams along rows
dist = sum(distTable(seqInds), 2);