function sequence = simSeq(uniqueNGrams, counts, seqLength, startState)

% SIMSEQ simulates a new sequence of length seqLength based on the input 
% matrix of unique n-grams and corresponding count vector, which define the
% empirical distribution that will be drawn on to generate the simulated 
% sequence.
% 
% Input
%   uniqueNGrams - an N x n matrix of n-grams observed in the set of
%                  training sequences (N is the number of unique n-grams)
%   counts       - an N x 1 vector of counts (the number of times each
%                  unique n-gram was observed in the training set)
%   seqLength    - the length of the simulated sequence
%   startState   - a 1 x (n - 1) vector defining the starting state
% 
% Output
%   sequence     - a 1 x seqLength sequence of states
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



% get n-gram order
n = size(uniqueNGrams, 2);

% check arguments
if size(uniqueNGrams, 1) ~= length(counts)
    error('Sizes of uniqueNGrams and counts do not match.')
end
[testGrams, ~] = countUniqueRows(uniqueNGrams);
if any(testGrams(:) ~= uniqueNGrams(:))
    error('uniqueNGrams must have only unique, sorted rows.')
end
if length(startState) ~= n - 1
    error('startState must be 1 element shorter than input n-grams.')
end

% get the largest state label (note, it doesn't matter if the number of
% possible states in the training data is larger than this in principle
% because we are only sampling from the states that were actually observed)
maxState = max(uniqueNGrams(:));

% get the unique (n - 1)-grams from uniqueNGrams
[roots, ~] = countUniqueRows(uniqueNGrams(:, 1:end - 1));

% fill in the transition probability matrix
transProb = zeros(size(roots, 1), maxState);
currentRoot = uniqueNGrams(1, 1:end-1);
jj = 1;
for ii = 1:size(uniqueNGrams, 1)
    % check to see if the root has changed
    if any(currentRoot ~= uniqueNGrams(ii, 1:end-1))
        % convert transProb from counts to probabilities
        transProb(jj, :) = transProb(jj, :) / sum(transProb(jj, :));
        
        % update current root
        currentRoot = uniqueNGrams(ii, 1:end-1);
        
        % increment transProb index
        jj = jj + 1;
    end
    
    % fill in the current row of the transition probability matrix
    transProb(jj, uniqueNGrams(ii, end)) = counts(ii);
end

% also normalise the last row
transProb(end, :) = transProb(end, :) / sum(transProb(end, :));

% get the cummulative probability for each row
cumProb = cumsum(transProb, 2);

% initialise the output sequence
sequence = NaN(1, seqLength);
sequence(1:n - 1) = startState;

% get random numbers for selecting next part of sequence
randNums = rand(seqLength, 1);

% generate the sequence
for ii = n:seqLength - length(startState)
    % find the current root
    [~, matchInd] = ismember(roots, sequence(ii-n+1:ii-1), 'rows');
    matchInd = logical(matchInd);
    
    % find the next element in the sequence
    [~, binInd] = histc(randNums(ii), cumProb(matchInd, :));
    
    % add 1 to binInd to change from zero indexing
    sequence(ii) = binInd + 1;
end