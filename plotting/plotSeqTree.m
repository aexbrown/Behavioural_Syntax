function plotSeqTree(compVec, gExp, ind)

% PLOTSEQTREE plots a tree structure indicating how the elements of gExp
% relate to each other.  compVec is also required so that the branches of
% the tree can be labelled with the number of occurances of the
% sub-sequences that we used to generate the tree when the original
% sequence was compressed.
%
% Input
%   compVec      - the compressed vector
%   gExp         - a number of rules by 1 cell array.  It is the expanded
%                  version of the input grammar obtained by applying the rules
%                  encoded in the grammar.  It contains the progressively
%                  expanded version of each of the rules in the grammar.  For
%                  any given entry gExp{ind}{end} contains only terminal
%                  symbols.
%   ind          - the index of the grammar in gExp that should be plotted

% loop through gExp
figure
xlim([0, numel(gExp{ind}{end}) + 1])
ylim([0, numel(gExp{ind}{end}) + 1.5])

% first plot the state labels of the terminal indices
for jj = 1:numel(gExp{ind}{end})
    text(jj, 0.5, num2str(gExp{ind}{end}(jj)), ...
        'HorizontalAlignment', 'center', 'FontSize', 20)
end

% loop to higher levels in gExp starting from the terminals
xCoords = 1:numel(gExp{ind}{end});
y = 0;
for jj = numel(gExp{ind}):-1:2
    % shift the y level
    y = y + 1;
    
    % plot vertical lines
    for kk = 1:numel(xCoords)
        line([xCoords(kk), xCoords(kk)], [y, y+1], 'LineWidth', 3)
    end
    
    % find the location of the new state following the merge and the
    % length difference of the sequence before and after merging
    newState = setdiff(gExp{ind}{jj-1}, gExp{ind}{jj});
    mergeStarts = ...
        find(gExp{ind}{jj-1} == newState);
    lengthDiff = (numel(gExp{ind}{jj}) - numel(gExp{ind}{jj-1})) / ...
       numel(mergeStarts);
    
    % draw the horizontal lines
    indShift = 0;
    for kk = 1:numel(mergeStarts)
        % every merged state shifts subsequent starts lengthDiff to the 
        % right.
        line([xCoords(mergeStarts(kk) + indShift), ...
            xCoords(mergeStarts(kk) + indShift + lengthDiff)], ...
            [y+1, y+1], 'LineWidth', 3)
        indShift = indShift + 1;
    end
    
    % update x coordinates
    for kk = 1:numel(mergeStarts)
        newX = xCoords(mergeStarts(kk)) + (xCoords(mergeStarts(kk) + lengthDiff) - xCoords(mergeStarts(kk)))/2;
        xCoords(mergeStarts(kk):mergeStarts(kk) + lengthDiff) = newX;
        xCoords = unique(xCoords);
    end
    
    % add text indicating the number of times the sequence
    % corresponding to the current branch was found during compression
    termInds = getTerminalInds(compVec, gExp, newState);
    newStateCounts = sum(termInds);
    
    for kk = 1:numel(mergeStarts)
        text(xCoords(kk)+0.1, y+1.2, num2str(newStateCounts), ...
            'FontSize', 20, 'Color', 'r', 'HorizontalAlignment', 'left')
    end
end

% add a final vertical line
xFinal = xCoords(1) + (xCoords(end) - xCoords(1))/2;
line([xFinal, xFinal], [y+1, y+2], 'LineWidth', 3)

% remove tick marks/axes
axis off
