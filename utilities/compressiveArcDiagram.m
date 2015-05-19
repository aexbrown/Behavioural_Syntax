function compressiveArcDiagram(compVec, gExp, sampRate, color, plotMode)

% COMPRESSIVEARCDIAGRAM plots an arc diagram that connects each repeated
% instance of a grammar rule found by the compressive sequence compression
% algorithm.  This is not the same rule used in the original arc diagram
% paper but the plot type is the same.
%
% Input
%   compVec  - the compressed vector
%   gExp     - a number of rules by 1 cell array.  It is the expanded
%              version of the input grammar obtained by applying the rules
%              encoded in the grammar.  It contains the progressively
%              expanded version of each of the rules in the grammar.  For
%              any given entry gExp{ii}{end} contains only terminal
%              symbols.
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



% loop through gExp from the highest level downwards
xMax = 1;
for ii = size(gExp, 1):-1:1
    disp(ii/size(gExp, 1))
    % get the terminal indices for the current rule
    termInds = getTerminalInds(compVec, gExp, ii + gExp{1}{1} - 1);
    termInds = find(termInds);
    
    % get the start and end indices for each group
%     ruleSize = size(gExp{ii - gExp{1}{1} + 1}{end}, 2);
    ruleSize = size(gExp{ii}{end}, 2);
    termStarts = termInds(1:ruleSize:end);
    termEnds = termInds(ruleSize:ruleSize:end);
    
    if max(termEnds) > xMax
        xMax = max(termEnds);
    end
    
    % loop through starts and plot the corresponding arc for each one
    for jj = 1:length(termStarts) - 1
        % plot the current arc
        plotArc([termStarts(jj), termEnds(jj), ...
            termStarts(jj+1), termEnds(jj+1)], color, ...
            sampRate, plotMode)
    end
end

% plot a baseline
line([1 xMax], [0 0], 'Color', [0 0 0], 'LineWidth', 3)

% make sure semicircles appear as semicircles
axis equal

% don't display axes
axis off

% add a scale bar that is 15% of the total length
line([xMax - xMax*0.15, xMax], [-xMax*0.05, -xMax*0.05], ...
    'Color', [0 0 0], 'LineWidth', 7)
text(xMax - xMax*0.15, -xMax*0.1, ...
    [num2str(round(xMax*0.15), 3) ' states'],...
    'FontSize', 14, 'FontWeight', 'bold')