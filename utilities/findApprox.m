function matchInds = findApprox(traj, patternTraj, r2Thresh)

% FINDAPPROX searches a matrix of trajectories for instances that
% approximately match the pattern trajectory.  Matches are defined as those
% with an R^2 value between the test trajectory and the pattern trajectory
% that is greater than or equal to r2Thresh, and where the sign of R is
% positive. (i.e. negative correlations shouldn't count as matches).
% 
% Input
%   traj        - A length(trajectory) by number of trajectories matrix
%   patternTraj - A length(trajectory) by 1 vector
%   r2Thresh    - The threshold R^2 value for definining a match
% 
% Output
%   matchInds   - The indices of any matching trajectories
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



% intialise
matchInds = zeros(1, size(traj, 2));

% loop through trajectories in traj
for ii = 1:size(traj, 2)
    % calculate the correlation between the current trajectory and the
    % pattern trajectory
    R = corrcoef(traj(:, ii), patternTraj);
    R2 = R(1, 2)^2;
    
    % check for a match
    if R2 >= r2Thresh && sign(R(1, 2)) > 0
        matchInds(ii) = 1;
    end
end