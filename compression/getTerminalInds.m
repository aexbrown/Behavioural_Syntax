function terminalInds = getTerminalInds(compVec, gExp, nonTerminal)

% GETTERMINALINDS calculates the indices of the terminals in the original
% uncompressed vector compVec was derived from corresponding to the given
% non-terminal.  You can't simply search for all instances of the lowest
% level of the expanded rule directly in the uncompressed vector because 
% some may not have been picked up in the original compression because of 
% the order in which it was carried out.
% 
% Input
%   compVec      - the compressed vector
%   gExp         - a number of rules by 1 cell array.  It is the expanded
%                  version of the input grammar obtained by applying the rules
%                  encoded in the grammar.  It contains the progressively
%                  expanded version of each of the rules in the grammar.  For
%                  any given entry gExp{ii}{end} contains only terminal
%                  symbols.
% 
% Output
%   terminalInds - A vector with the start indices of the terminals
%                  **in the uncompressed vector** corresponding to the
%                  input non-terminal
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



% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% this funciton uses strrep on a vector of doubles, which works but would
% normally issue a warning saying it doesn't work...  It's called many
% times, so turn off warning
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
warning('off','MATLAB:strrep:InvalidInputType')


parentRepVec = compVec;
while any(parentRepVec >= gExp{1}{1})
    % get the current maximum non-terminal
    maxNT = max(parentRepVec);
    
    % if maxNT is the non-terminal of interest, do the replacements with
    % negative copies of itself (the negative is just to ensure the values
    % are less than gExp{1}{1} so that the while loop continues)
    if maxNT == nonTerminal
        parentRepVec = strrep(parentRepVec, maxNT, ...
            repmat(-maxNT, 1, numel(gExp{maxNT - gExp{1}{1} + 1}{end})));
    else
        % this is not the non-terminal we're searching for.  Expand it by
        % one level and continue search.
        parentRepVec = strrep(parentRepVec, maxNT, ...
            gExp{maxNT - gExp{1}{1} + 1}{2});
    end
end

% find the indices of the given non-terminal
terminalInds = parentRepVec == -nonTerminal;

% turn warning back on
warning('on','MATLAB:strrep:InvalidInputType')