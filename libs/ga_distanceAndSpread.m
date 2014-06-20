function [avgDistance, spread] = ga_distanceAndSpread(Distance,Rank,score,score_old)
%distanceAndSpread Calculate average distance and spread of a population
%   distanceAndSpread is private to GAMULTIOBJ

%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2009/08/29 08:25:54 $


% Get individual and their distance measure from first front
distance = Distance((Rank == 1));
% Remove infinite distances
distance = distance(isfinite(distance));
if ~isempty(distance)
    num_solution = length(distance);
    avg_distance = sum(distance)/num_solution;
    avgDistance = norm(distance - avg_distance)/sqrt(num_solution);
else
    avg_distance = 0;
    avgDistance = 0; % All solutions are at extreme
end

[unused,index] = min(score,[],1);
extremeParetoSol_new = score(index,:);

[unused,index] = min(score_old,[],1);
extremeParetoSol_old = score_old(index,:);

% Calculate average distance between extreme solution
extremeParetoDistance = 0;
for i = 1:size(extremeParetoSol_new,1)
    extremeParetoDistance = extremeParetoDistance + norm((extremeParetoSol_old(i,:) - extremeParetoSol_new(i,:)));
end

if ~isfinite(extremeParetoDistance)
   extremeParetoDistance = 0; 
end
% Calculate spread of the non-dominated front
if size(extremeParetoSol_new,2) > 1 && (extremeParetoDistance > 0 ...
        || avgDistance > 0)
    spread = (extremeParetoDistance + avgDistance)/ ...
        (extremeParetoDistance + size(extremeParetoSol_new,2)*avg_distance);
else
    spread =  extremeParetoDistance;
end
