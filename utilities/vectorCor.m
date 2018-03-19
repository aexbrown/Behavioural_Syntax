function Rmat = vectorCor(mat1, mat2)

% VECTOR calculates the correlation coefficient between the rows of mat1
% and mat2 and returns the result as a matrix (analogous to a distance
% matrix).
% 
% Input
%   mat1,2 - Two matrices with the same number of columns.
% 
% Output
%   Rmat   - A matrix of correlation coefficients.  Rmat(i, j) is the
%            correlation coefficient between the ith row of mat1 and the
%            jth row of mat2.  Rmat has dimensions size(mat1, 1) by
%            size(mat2, 1).
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




Rmat = NaN(size(mat1, 1), size(mat2, 1));
for ii = 1:size(mat1, 1)
    for jj = 1:size(mat2, 1)
        % calculate the correlation
        R = corrcoef(mat1(ii, :), mat2(jj, :));
        Rmat(ii, jj) = R(1, 2);
    end
end