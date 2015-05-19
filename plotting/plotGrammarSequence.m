function plotGrammarSequence(gExp, postures, sequenceInds)

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
for ii = 1:length(sequenceInds)
    currentInd = sequenceInds(ii);
    currentSeq = gExp{currentInd}{end};
    figure
    hold on
    for jj = 1:length(currentSeq)
        % convert the current matching posture to xy-coordinates
        [x, y] = angle2skel(postures(currentSeq(jj), :)', ...
            -pi/2, arclength);
        
        % plot the matching posture over the original skeleton
        plot(x - mean(x) + jj * offsetFactor, ...
            y - mean(y), 'Color', 'r', 'LineWidth', 3)
        
        % mark the head
        plot(x(1) - mean(x) + jj * offsetFactor, ...
            y(1) - mean(y), 'Color', 'b', 'LineWidth', 5)
    end
    axis equal
    hold off
end