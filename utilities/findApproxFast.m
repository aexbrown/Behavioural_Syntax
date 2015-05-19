function matchInds = findApproxFast(nGrams, pattern, r2Thresh, rMat)

% FINDAPPROXFAST searches a matrix of trajectories for instances that
% approximately match the pattern trajectory.  Matches are defined as those
% with an R^2 value between each part of the test trajectory and each part
% of the the pattern trajectory that is greater than or equal to r2Thresh, 
% and where the sign of R is positive. (i.e. negative correlations 
% shouldn't count as matches).
% 
% SEE ALSO: findApprox.m
% 
% Input
%   nGrams      - an N by n matrix of n-grams (each row is an n-gram)
%   pattern     - a 1 by n vector to match to the rows of nGrams
%   r2Thresh    - The threshold R^2 value for definining a match
%   rMat        - A numPostures by numPostures matrix that includes the R
%                 values for each posture against each other posture.
% 
% Output
%   matchInds   - The indices of any matching trajectories
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



% calculate the R^2 and sign matrices
r2Mat = rMat.^2;
signMat = sign(rMat);

% convert nGrams to indices for getting elements from r2Mat
inds = sub2ind(size(r2Mat), nGrams, repmat(pattern, size(nGrams, 1), 1));
r2Comparisons = r2Mat(inds);
signComparisons = signMat(inds);

% check the R2 values and the signs
combinedCheck = r2Comparisons >= r2Thresh & signComparisons > 0;
matchInds = all(combinedCheck');

% % ----debugging, check result with a loop-----
% % intialise
% matchInds2 = zeros(1, size(nGrams, 1));
% 
% % loop through comparison values
% for ii = 1:size(nGrams, 1)
%     % check the signs
%     if all(signComparisons(ii, :) > 0)
%         % check the fit
%         if all(r2Comparisons(ii, :) >= r2Thresh)
%             matchInds2(ii) = 1;
%         end
%     end
% end
% 
% if any(matchInds ~= matchInds2)
%     disp('Something is wrong...')
% end