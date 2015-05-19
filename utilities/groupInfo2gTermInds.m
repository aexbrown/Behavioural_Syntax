function gTermInds = groupInfo2gTermInds(groupInfo, postureNum, gTerminals)

% convert a list of column names from a clustergram group information
% struct, to the corresponding locations in a gTerminals cell.  This of
% course only works if the gTerminals cell was used to create the
% clustergram object where the groupInfo comes from.
% 
% Input
%   groupInfo - a struct that can be saved to the workspace manually by
%               right-clicking on nodes in a clustergram
% 
% Output
%   gTermInds - the indices in the gTerminals cell corresponding to the
%               column names in the clustergram
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



% manually select groups from the clustergram and export group information
% to the workspace
groupInds = str2num(char(groupInfo.ColumnNodeNames));

% convert group indices to corresponding entries in gTerminals.
% first get any direct posture indices
gTermInds = groupInds(groupInds <= postureNum);

% get the rest of the indices.  Use mod to reference back to elements in
% gTerminals
gTermInds = [gTermInds; ...
    mod(groupInds(groupInds > postureNum), numel(gTerminals))];