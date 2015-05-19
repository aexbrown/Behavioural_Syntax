function accCurve = nGramAccumulateNumerical(nGrams)

% NGRAMACCUMULATENUMERICAL Works sequentially through the input matrix of 
% n-grams and determines where unique occurrances arise.  These are output
% in an accumulation curve.
%
% Input
%   nGrams   - A matrix of n-grams
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



% get the unique nGrams
uGrams = unique(nGrams, 'rows');

% set flags for each of the unique n-grams
flags = zeros(1, size(uGrams, 1));

% set flag for each of the input n-grams
isnew = zeros(1, size(nGrams, 1));

for ii = 1:size(nGrams, 1)
    % get the index of the current n-gram in the unique list
    ind = all(uGrams == repmat(nGrams(ii, :), size(uGrams, 1), 1), 2);
    
    % has the current n-gram been seen before?
    if ~flags(ind)
        % if not, update flags
        isnew(ii) = 1;
        flags(ind) = 1;
    end
end

% calculate the accumulation
accCurve = cumsum(isnew);