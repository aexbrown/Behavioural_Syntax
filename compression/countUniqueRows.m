function [uniqueRows, counts] = countUniqueRows(data)

% COUNTUNIQUEROWS counts the number of occurances of each of the unique
% rows found in dataMat.
% 
% Input
%   data       - A matrix or cell array of row vectors whose unique rows 
%                will be counted
% 
% Output
%   uniqueRows - A matrix containing only the unique rows of dataMat
%   counts     - A vector of counts.  The first element of counts gives the
%                number of occurances in dataMat of the corresponding row
%                in uniqueRows
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



% check whether the input is a cell or matrix
if iscell(data)
    % convert the cell to a cell of strings for use with unique
    dataString = cellfun(@num2str, data, 'UniformOutput', false);
    
    % get the unique rows and their counts
    [uniqueRows, ~, c] = unique(dataString);
    counts = hist(c, size(uniqueRows, 1));
    
    % convert uniqueRows back to a cell of row vectors
    uniqueRows = cellfun(@str2num, uniqueRows, 'UniformOutput', false);
else
    % get the unique rows and their counts
    [uniqueRows, ~, c] = unique(data, 'rows');
    counts = hist(c, size(uniqueRows, 1));
end