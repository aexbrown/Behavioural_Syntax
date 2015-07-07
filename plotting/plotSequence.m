function plotSequence(postures, sequence, varargin)

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




if nargin == 2
    bodyColor = 'r';
    headColor = 'b';
    width = 2;
    offsetX = 0;
    offsetY = 0;
elseif nargin == 3
    bodyColor = varargin{1};
    headColor = 'b';
    width = 2;
    offsetX = 0;
    offsetY = 0;
elseif nargin == 4
    bodyColor = varargin{1};
    headColor = varargin{2};
    width = 2;
    offsetX = 0;
    offsetY = 0;
elseif nargin == 5
    bodyColor = varargin{1};
    headColor = varargin{2};
    width = varargin{3};
    offsetX = 0;
    offsetY = 0;
elseif nargin == 6
    bodyColor = varargin{1};
    headColor = varargin{2};
    width = varargin{3};
    offsetX = varargin{4};
    offsetY = 0;
elseif nargin == 7
    bodyColor = varargin{1};
    headColor = varargin{2};
    width = varargin{3};
    offsetX = varargin{4};
    offsetY = varargin{5};
end

arclength = 1;
offsetFactor = 0.3;

hold on
for jj = 1:length(sequence)
    % convert the current matching posture to xy-coordinates
    [x, y] = angle2skel(postures(sequence(jj), :)', ...
        -pi/2, arclength);
    
    % plot the matching posture over the original skeleton
    line(x - mean(x) + jj * offsetFactor + offsetX, ...
        y - mean(y) + offsetY, 'Color', bodyColor, 'LineWidth', width)
    
    % mark the head
    plot(x(1) - mean(x) + jj * offsetFactor + offsetX, ...
        y(1) - mean(y) + offsetY, ...
        '.', 'MarkerSize', width*15, 'Color', headColor)
end
ylim([-0.6, 0.6])
axis equal
hold off