function gTerminals = getGrammarTerminals(grammar)

% GETGRAMMARTERMINALS returns the only the terminal symbols of the expanded
% version of the input grammar.
%
% Input
%   grammar    - a number of rules by 2 cell array. The first column has
%                the left hand side of each replacement rule while the
%                second column has the right hand side (so the first column
%                lists all non-terminals in the grammar).
%
% Output
%   gTerminals - a number of rules by 1 cell array containing the expanded
%                terminal symbols of the input grammar.
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



% expand the grammar
[gExp, ~] = expandGrammar(grammar, []);

% extract the terminal symbols
gTerminals = cell(size(gExp));
for ii = 1:size(gExp, 1)
    % add only the terminals to gTerminals
    gTerminals{ii} = gExp{ii}{end};
end