function boo15()
global STACK;
% Only append MSE module if needed
mods = find_modules(STACK, 'mean_squared_error', true);
if isempty(mods)
    nmse();    
end

% Very rapid step growth, but absolute
fit_boo(struct('StopAtSeconds', 240, ...
               'StopAtStepNumber', 10000, ...
               'StopAtStepSize', 10^-12, ...
               'StopAtRelScoreDelta', 10^-12, ...
               'StopAtAbsScoreDelta', 10^-12, ...
               'StepSize', 1.0, ...
               'StepRel', false, ...
               'StepRelRecalcEvery', 1, ...
               'StepRelMin', 10^-3, ...
               'StepRelMax', 10^3, ...
               'StepGrowth', 2, ...
               'StepShrink', 0.5, ...
               'Elitism', false, ...
               'EliteParams', 5, ...
               'EliteSteps', 5));