function plotSequenceDensity(postures, sequence1, sequence2, offset)

% varargin are headColor, bodyColor, and width. Colors can be either one of
% the allowed strings (e.g. 'r' or 'b') or a 1x3 RGB vector.  Width is a
% number specifying the line width
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



arclength = 1;
offsetFactor = 0.3;

% prepare the first sequence points
xTotal = [];
yTotal = [];
for ii = size(sequence1, 1)
    for jj = 1:size(sequence1, 2)
        % convert the current matching posture to xy-coordinates
        [x, y] = angle2skel(postures(sequence1(jj), :)', ...
            -pi/2, arclength);
        
        
        % offset coordinates
        x = x - mean(x) + jj * offsetFactor;
        y = y - mean(y);
        
        xTotal = [xTotal, x];
        yTotal = [yTotal, y];
    end
end

% prepare the second sequence points
for ii = size(sequence2, 1)
    for jj = 1:size(sequence2, 2)
        % convert the current matching posture to xy-coordinates
        [x, y] = angle2skel(postures(sequence2(jj), :)', ...
            -pi/2, arclength);
        
        
        % offset coordinates
        x = x - mean(x) + jj * offsetFactor + offset;
        y = y - mean(y);
        
        xTotal = [xTotal, x];
        yTotal = [yTotal, y];
    end
end

ndhist(xTotal, yTotal, 'max', 'bins', 15, 'filter', 13);
axis equal
