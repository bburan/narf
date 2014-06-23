function state = ga_stepgamultiobj(subpopIndex,thisPopulation,options,state,GenomeLength,FitnessFcn)
%STEPGAMULTIOBJ perform one step using a variant of NSGA-II algorithm
%   This function is private to GAMULTIOBJ.

%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2009/08/29 08:27:21 $

popSize    = size(state.Population(thisPopulation,:),1);
population = state.Population(thisPopulation,:);
score      = state.Score(thisPopulation,:); 
rank       = state.Rank(thisPopulation,:);
Distance   = state.Distance(thisPopulation,:);
score_old  = score;
numObj = size(score,2);

% How many crossover offspring will there be from each source?
nXoverKids = round(options.CrossoverFraction * popSize);
nMutateKids = popSize - nXoverKids;
% how many parents will we need to complete the population?
nParents = 2 * nXoverKids + nMutateKids;


% Selection.
parents = feval(options.SelectionFcn,[-rank,Distance],nParents,options,options.SelectionFcnArgs{:});
% Shuffle to prevent locality effects. 
parents = parents(randperm(length(parents)));

% Everyones parents are stored here for genealogy display
state.Selection = [parents'];

% Here we make all of the members of the next generation
xoverKids  = feval(options.CrossoverFcn,parents(1:(2 * nXoverKids)),options,GenomeLength, ...
    FitnessFcn,score,population,options.CrossoverFcnArgs{:});
mutateKids = feval(options.MutationFcn,parents((1 + 2 * nXoverKids):end),options,GenomeLength, ...
    FitnessFcn,state,score,population,options.MutationFcnArgs{:});

% Group them into the next generation
nextPopulation = [xoverKids ; mutateKids ];

% Score the population
if strcmpi(options.Vectorized, 'off') 
    nextScore = fcnvectorizer(nextPopulation,FitnessFcn,numObj,options.SerialUserFcn);
else
    [nextScore] = FitnessFcn(nextPopulation);
end

% Update the current population
% state.Population(thisPopulation,:) = nextPopulation; % this is idiotic
% state.Score(thisPopulation,:) = nextScore; % same
% Update function evaluation counter
state.FunEval = state.FunEval + size(nextScore,1);

%--Prepare for next generation--

% Combine new and old population
population = [population;nextPopulation];
score = [score; nextScore];

% addition: will update the diversity score
n_sol = size(population,1);
dists = squareform(pdist(population));
score(:,2) = -min(dists+max(max(dists))*eye(n_sol));
% we remove the solutions that are really too bad:
toremove = find(score(:,1) > 5 * min(score(:,1)), 1);
if ~isempty(toremove)
    if length(toremove) >= popSize
        toremove = toremove(1:(length(toremove)-popSize-1));
    end
    score = score(~ismember(1:(2*popSize), toremove), :);
    population = population(~ismember(1:(2*popSize), toremove), :);
end

if size(score,1)< popSize
    fprintf('ouille.\n')
    fprintf('current pop size = %d (%d)\n',size(score,1),size(population,1));
end

% fprintf('best indiv score = %f\n',min(score(:,1)));

% Sort combined population and pick best 'popSize' individuals for next
% generation
[state.Population(thisPopulation,:),state.Score(thisPopulation,:), ...
    state.Rank(thisPopulation,:),state.Distance(thisPopulation,:)] = ...
    ga_rankAndDistance(population,score,options,popSize);

% Calculate average distance and spread for the Pareto front
[state.AverageDistance(subpopIndex), state.Spread(state.Generation+1,subpopIndex)] = ...
    ga_distanceAndSpread(state.Distance(thisPopulation,:),state.Rank(thisPopulation,:), ...
    state.Score(thisPopulation,:),score_old);

