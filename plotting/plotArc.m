function plotArc(X, color, sampRate, plotMode)

% PLOTARC plots a semicircular arc between two line segments defined by the
% four points in X.
% 
% See SEMICIRCLE, defined below
% 
% Input
%   X        - A 1x4 vector that defines the arc edges.  The arc is a patch
%              defined by a semi-circle connecting X(1) -> X(4), a straight
%              line from X(4) -> X(3), and a semi-circle from X(3) -> X(2).
%   color    - A 1x4 vector defining the color and transparency of the
%              patch in the following format: [r g b alpha]
%   sampRate - The sampling rate for creating the arc.  If sample rate is
%              1, there will be the same number of points as the arc
%              diameter.  Setting sampRate less than 1 will give
%              proportionately fewer points.
%   plotMode - 'patch' or 'line'.  'patch' only works well if width of the
%              arc is not too much smaller than the radius.  Otherwise it
%              just looks like a line anyway.  In this case 'line' is
%              better because it uses fewer points and the matlab plotter
%              renders it better on screen.
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



if strcmp(plotMode, 'patch')
% get the patch coordinates
[firstArcX, firstArcY] = semicircle(X(1), X(4), sampRate);
[secondArcX, secondArcY] = semicircle(X(2), X(3), sampRate);

% plot the patch
patch([firstArcX, secondArcX(end:-1:1)], ...
    [firstArcY, secondArcY(end:-1:1)], color(1:3), ...
    'FaceAlpha', color(4), 'LineStyle', 'none')
elseif strcmp(plotMode, 'line')
    % get line coordinates
    [arcX, arcY] = ...
        semicircle(X(1) + (X(2) - X(1))/2, X(3) + (X(4) - X(3))/2,...
        sampRate);

    % plot the line with a width proportional to the difference between
    % X(1) and X(2)
    line(arcX, arcY, 'LineWidth', (X(2) - X(1))*0.5, 'Color', color(1:3))
%     patch([arcX NaN], [arcY NaN], color(1:3), ...
%         'FaceAlpha', color(4), 'LineStyle', 'none')
else
    error('plotMode must be either patch or line')
end


function [x, y] = semicircle(startX, endX, sampRate)

% SEMICIRCLE generates points on a semicircle between startX and endX with
% a number of points equal to the diameter*sampRate (the sample rate).

radius = (endX - startX)/2;
xCenter = startX + radius;
theta = 0 : pi/(sampRate*2*radius) : pi;
x = radius * cos(theta) + xCenter;
y = radius * sin(theta);