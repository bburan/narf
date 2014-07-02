function [pop,score,nonDomRank,Distance]  = ga_rankAndDistance(pop,score,options,nParents)
%rankAndDistance Assign rank and distance measure to each individual

%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2009/08/29 08:27:00 $

   
if nargin < 4
    nParents = size(pop,1);
end

ParetoFraction = options.ParetoFraction;
nScore = size(score,2);
if nScore == 1 % single objective
    nonDomRank = ga_nonDominatedRank(score,nParents);
    % Remove individuals with infinite rank
    index = isinf(nonDomRank);
    nonDomRank(index) = [];
    pop(index,:) = [];
    score(index,:) = [];
else
    nonDomRank = ga_nonDominatedRank(score);
end
popSize = size(pop,1);
Distance = zeros(popSize,1);
numRank = unique(nonDomRank); % numRank will be sorted

% Compute crowding distance for individuals in each front
for i = numRank'
   % Get individual from each front
   index = (nonDomRank == i);
   Distance(index) = options.DistanceMeasureFcn(pop(index,:),score(index,:),options,options.DistanceMeasureFcnArgs{:}); 
end

% If populations were not combined then no need to trim the population
if nParents == popSize
    % do nothing
else
    [pop,score,nonDomRank,Distance] = ga_trimPopulation(pop,score,nonDomRank,Distance, ...
        popSize,nScore,nParents,ParetoFraction);
end