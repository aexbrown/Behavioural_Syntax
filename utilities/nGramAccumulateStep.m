function accCurve = nGramAccumulateStep(nGrams, stepSize)

% NGRAMACCUMULATENUMERICAL Works through the input matrix of n-grams in
% chunks to determine growth of number of unique sequences with growing
% number of observed n-grams.  It can be much faster that 
% nGramAccumulateNumerical if the step size is large.
%
% Input
%   nGrams   - A matrix of n-grams
%   stepSize - The step size for moving through n-grams.
%
% Output
%   accCurve - The accumulation curve of previously unseen n-grams as they
%              occur in nGrams scanned from start to finish
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



% get the chunk boundaries
chunkBounds = stepSize:stepSize:size(nGrams, 1);

% initialise
accCurve = [1, NaN(1, numel(chunkBounds))];

% loop through n-gram chunks
for ii = 1:numel(chunkBounds)
    accCurve(ii+1) = size(unique(nGrams(1:chunkBounds(ii), :), 'rows'), 1);
end