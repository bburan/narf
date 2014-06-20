function [pop,score,nonDomRank,crowdingDistance] = trimPopulation(pop,score,nonDomRank,crowdingDistance, ...
                                                   popSize,nScore,nParents,ParetoFraction)
%trimPopulation Trim population to get nParents with specified distribution
%among all fronts

%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2009/08/29 08:27:26 $


if nScore > 1 % front diversity make sense for multi-objective
    [pop,score,nonDomRank,crowdingDistance,popSize] = frontDiversity(pop,score,nonDomRank,crowdingDistance, ...
        popSize,nParents,ParetoFraction);
end
% trim individuals to nParents size based on distance
[pop,score,nonDomRank,crowdingDistance] = distanceDiversity(pop,score,nonDomRank,crowdingDistance, ...
    nParents,popSize);

%-------------------------------------------------------------
function [pop,score,nonDomRank,crowdingDistance] = distanceDiversity(pop,score,nonDomRank, ...
    crowdingDistance,nParents,popSize)
% Remove indivduals from last front if length(nonDomRank) is greater than nParents
to_remove = length(nonDomRank) - nParents;
if to_remove > 0
    [unused,I] = sortrows([nonDomRank, crowdingDistance],[1,-2]);
    indiv_to_remove = I(popSize:-1:(popSize - to_remove + 1));
    pop(indiv_to_remove,:) = [];
    score(indiv_to_remove,:) = [];
    nonDomRank(indiv_to_remove) = [];
    crowdingDistance(indiv_to_remove) = [];
end

%-------------------------------------------------------------------------------
function [pop,score,nonDomRank,crowdingDistance,popSize] = frontDiversity(pop,score,nonDomRank,crowdingDistance, ...
                                                   popSize,nParents,ParetoFraction)
numRank = unique(nonDomRank); % numRank will be sorted
% Trim front size based on a pre-determined count
totalNumOfRank = length(numRank);
retainIndiv = zeros(totalNumOfRank,1); % Individuals to be retained in each front
availableIndiv = zeros(totalNumOfRank,1); % Available individuals in each front

% First rank is treated differently for number of desired individuals
availableIndiv(1) = nnz(nonDomRank == 1);
allowedIndiv = ceil(nParents*ParetoFraction);
if allowedIndiv <= availableIndiv(1)
    retainIndiv(1) = allowedIndiv;
else
    retainIndiv(1) =  availableIndiv(1);
end

% If there are more than one fronts then we calculate the maximum number of
% individuals we can retain
if totalNumOfRank > 1
    gpRatio = 0.8; % Geometric progression ratio
    carryover = 0;
    for i = numRank(2:end)'
        availableIndiv(i) = nnz(nonDomRank == i);
        % allowedIndiv redcuces as geometric progression
        allowedIndiv = ceil(nParents*gpRatio^(i-1)*((1-gpRatio)/(1-gpRatio^totalNumOfRank))) + carryover;
        if allowedIndiv <= availableIndiv(i)
            retainIndiv(i) = allowedIndiv;
        else
            retainIndiv(i) =  availableIndiv(i);
            carryover = allowedIndiv - availableIndiv(i);
        end
    end
end

% Adjust retainIndiv if required
front = totalNumOfRank; increase = 0.10;
while sum(retainIndiv) <= ceil(nParents)
    if retainIndiv(front) < availableIndiv(front)
        retainIndiv(front) = min(availableIndiv(front),ceil(retainIndiv(front)*(1+increase)));
        increase = increase*0.95;
    end
    if front == 1
        front = totalNumOfRank;
        increase = 0.10;
    else
        front = front - 1;
    end
end

for i = numRank'
   popIndices = 1:popSize;
   index = (nonDomRank == i);
   if retainIndiv(i) < nnz(index)   % We need to trim the front
       if retainIndiv(i) == 0
           continue;
       end
       % Select individuals after a tournament selection
       % Choose the players for each round
       playerlist = popIndices(index);
       losers = tournament(playerlist,crowdingDistance,retainIndiv(i)); 
       pop(losers,:) = [];
       score(losers,:) = [];
       nonDomRank(losers) = [];
       crowdingDistance(losers) = [];
       popSize = size(score,1);
   else
       % No need to modify the front; do nothing
   end
end

%-------------------------------------------------------------
function losers = tournament(playerlist,expectation,allowed)
%tournament between players based on their distance values

[sortedExpectation,sortedIndex] = sort(expectation(playerlist),'descend');
losers = playerlist(sortedIndex(allowed+1:end));
