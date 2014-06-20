function options = initialPopulationCheck(options)
%initialPopulationCheck Validates initial population

%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2009/10/10 20:09:21 $

% Make a *guess* at numberOfObjectives 
if ~isempty(options.InitialScores)
    numberOfObjectives = size(options.InitialScores,2);
else
    numberOfObjectives = 1; % safe for testing
end

% Scores for single objective should be a column vector
if numberOfObjectives == 1
   options.InitialScores = options.InitialScores(:);
end
size_initPop  = size(options.InitialPopulation,1);
[len_initScore,numObj] = size(options.InitialScores);

% No tests if initial pop and scores are empty
if size_initPop == 0 && len_initScore == 0
    return;
end

popSize  = sum(options.PopulationSize);

if size_initPop > popSize
    warning('globaloptim:initialPopulationCheck:initPopLength','Initial population contains more individuals than population size;ignoring extra individuals.');
    options.InitialPopulation(popSize+1:size_initPop,:) = [];
    size_initPop  = size(options.InitialPopulation);
end

if len_initScore > popSize
    warning('globaloptim:initialPopulationCheck:initScoreLength','Number of initial scores more than population size;ignoring extra scores.');
    options.InitialScores(popSize+1:len_initScore,:) = [];
    len_initScore = size(options.InitialScores,1);
end

if len_initScore > size_initPop
    warning('globaloptim:initialPopulationCheck:initScoreAndPopLength','Number of initial scores exceeds the number of individuals in the initial population.');
end

