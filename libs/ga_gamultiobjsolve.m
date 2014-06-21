function [x,fval,exitFlag,output,population,scores] = gamultiobjsolve(FitnessFcn,GenomeLength, ...
     Aineq,bineq,Aeq,beq,lb,ub,options,output)
%GAMULTIOBJSOLVE Genetic algorithm multi-objective solver.

%   Copyright 2007-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2009/10/10 20:09:17 $



% Create initial state: population, scores, status data
state = ga_gamultiobjMakeState(GenomeLength,FitnessFcn,output.problemtype,options);

currentState = 'init';
% Give the plot/output Fcns a chance to do any initialization they need.
% state = gadsplot(options,state,currentState,'Genetic Algorithm'); %%%%%%%
% [state,options] = gaoutput(FitnessFcn,options,state,currentState); %%%%%%

% Setup display header 
if  options.Verbosity > 1
    fprintf('\n                           Average            Average\n');
    fprintf('Generation   f-count    Pareto distance    Pareto spread\n');
end

historyBest = repmat(Inf, options.StallGenLimit, 1);
historyPointer = 1;

currentState = 'iter';
% Run the main loop until some termination condition becomes true
exitFlag = [];
while true
       state.Generation = state.Generation + 1;
        % check to see if any stopping criteria have been met
%        [state,exitFlag,reasonToStop] = ga_gamultiobjConverged(options,state);
       exitFlag = [];
       if state.Generation > options.Generations
           exitFlag = 1;
           reasonToStop = 'Maximal number of iterations reached => exiting\n';
       end
       if range(historyBest) < options.TolFun
           exitFlag = 1;         
           reasonToStop = 'Optimization stalled => exiting\n';
       end
       if ~isempty(exitFlag)
           break;
       end
       
       
            % Print a newline and important info every 10 iterations
            if isequal(mod(state.Generation-1, 5), 1)
                fprintf('\nGeneration %5d => current scores: [%f...%f...%f]\n', state.Generation, ...
                    min(state.Score(:,1)), mean(state.Score(:,1)), max(state.Score(:,1)));
                [~, best] = min(state.Score(:,1));
                fprintf('best indiv last parameter is %f\n',state.Population(best,end));
                % dbtickqueue(n_iters);
            end
       
        % Repeat for each sub-population (element of the PopulationSize vector)
        offset = 0;
        totalPop = options.PopulationSize;
        % Each sub-population loop
        for pop = 1:length(totalPop)
            populationSize =  totalPop(pop);
            thisPopulation = 1 + (offset:(offset + populationSize - 1));
            % Empty population is also possible
            if isempty(thisPopulation)
                continue; 
            end
            state = ga_stepgamultiobj(pop,thisPopulation,options,state,GenomeLength,FitnessFcn);
            offset = offset + populationSize;
        end
        
        historyBest(historyPointer) = min(state.Score(:,1));
        historyPointer = historyPointer + 1;
        if historyPointer > options.StallGenLimit
            historyPointer = 1;
        end
        
        % Migration
%         state = migrate(FitnessFcn,GenomeLength,options,state); %%%%%%%%%

        % Output and plot functions
%         state = gadsplot(options,%state,currentState,'Genetic Algorithm');
%         [state,options] = gaoutput(FitnessFcn,options,state,currentState);
end % End while loop
% Update output structure
output.generations = state.Generation;
output.message = reasonToStop;
fprintf(output.message);

% If sub-population model is used, merge all sub-population and perform
% another non-dominated sorting
if length(options.PopulationSize) > 1
    [state.Population,state.Score,state.Rank,state.Distance]  = ...
        rankAndDistance(state.Population,state.Score,options);
    % Calculate average distance and spread
    [output.averagedistance,output.spread] = distanceAndSpread(state.Distance, ...
        state.Rank,state.Score,state.Score);
else
    % Calculate front statistics for output structure
    output.averagedistance = state.AverageDistance;
    output.spread = state.Spread(end);
end

% Find and return the solutions on Pareto front
fval = state.Score(state.Rank == 1,:);
x = state.Population((state.Rank == 1),:);


% % A hybrid scheme; try another minimization method
% if ~isempty(options.HybridFcn)
%     if strcmpi(options.PopulationType,'doubleVector')
%         state  = gamultiobjHybrid(FitnessFcn,x,fval,state,Aineq,bineq,Aeq,beq,lb,ub,options);
%         % Calculate front statistics for output structure
%         [output.averagedistance,output.spread] = distanceAndSpread(state.Distance,state.Rank,state.Score,state.Score);
%         % Find and return the solutions on Pareto front
%         fval = state.Score(state.Rank == 1,:);
%         x = state.Population((state.Rank == 1),:);
%     else
%         warning('globaloptim:gamultiobjsolve:notValidHybrid','''HybridFcn'' can only be used with ''doubleVector'' population; ignoring ''HybridFcn'' option');
%     end
% end

output.funccount   = state.FunEval;
population = state.Population;
scores = state.Score;

currentState = 'done';
% Give the Output functions a chance to finish up
% gadsplot(options,state,currentState,'Genetic Algorithm');
% gaoutput(FitnessFcn,options,state,currentState);

