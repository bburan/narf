function boo09()
mse();

fit_boo(struct('StopAtSeconds', 120, ...
               'StopAtStepNumber', 10000, ...
               'StopAtStepSize', 10^-12, ...
               'StopAtRelScoreDelta', 10^-12, ...
               'StopAtAbsScoreDelta', 10^-12, ...
               'StepSize', 1.0, ...
               'StepRel', true, ...
               'StepRelRecalcEvery', 10, ...
               'StepRelMin', 10^-3, ...
               'StepRelMax', 10^3, ...
               'StepGrowth', 1.0, ...
               'StepShrink', 0.5, ...
               'Elitism', true, ...
               'EliteParams', 5, ...
               'EliteSteps', 10));