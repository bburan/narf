function boo05()
mse();
                     
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
               'StepGrowth', 1.4, ...
               'StepShrink', 0.5, ...
               'Elitism', true, ...
               'EliteParams', 5, ...
               'EliteSteps', 5));