function [gExp, vecExp] = expandGrammar(grammar, compVec)

% EXPANDEDGRAMMAR expands the nested structure of a grammar based on the
% rules it encodes.  Also creates a matrix that contains the progressively
% uncompressed versions of the compressed vector compVec.
%
% Input
%   grammar  - a number of rules by 2 cell array. The first column has the
%              left hand side of each replacement rule while the second
%              column has the right hand side (so the first column lists
%              all non-terminals in the grammar).
%   compVec  - a vector that has been compressed using grammar.
%
% Output
%   gExp     - a number of rules by 1 cell array.  It is the expanded
%              version of the input grammar obtained by applying the rules
%              encoded in the grammar.  It contains the progressively
%              expanded version of each of the rules in the grammar.  For
%              any given entry gExp{ii}{end} contains only terminal
%              symbols.
%   vecExp   - 
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



% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% this funciton uses strrep on a vector of doubles, which works but would
% normally issue a warning saying it doesn't work...  It's called many
% times, so turn off warning
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
warning('off','MATLAB:strrep:InvalidInputType')

% initialise expanded grammar and expanded vector
gExp = cell(size(grammar, 1), 1);
vecExp = cell(size(grammar, 1)+1, 1);
vecExp{size(grammar, 1)+1} = compVec;

% get the non-terminal of the first rule
smallestNT = grammar{1, 1};

% count down through grammar from the last entry and apply compression 
% rules in reverse
for ii = size(grammar, 1):-1:1
    % initialise the current entry of the expanded grammar
    gExp{ii} = grammar(ii, 1);
    
    % get the first non-terminal symbol to replace
    currentNT = grammar{ii, 1};
    
    % make the replacements in compVec
    vecExp{ii, 1} = strrep(vecExp{ii + 1, :}, currentNT, grammar{ii, 2});
    
    % make the replacements in the grammar until only terminals are left
    while any(gExp{ii}{end} >= smallestNT)
        % get the maximum current entry
        maxNT = max(gExp{ii}{end});
        
        % replace the maximum entry
        gExp{ii}{end+1} = strrep(gExp{ii}{end}, maxNT, ...
            grammar{maxNT - smallestNT + 1, 2});
    end
end

% turn warning back on
warning('on','MATLAB:strrep:InvalidInputType')