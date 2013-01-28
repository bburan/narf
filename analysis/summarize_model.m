function s = summarize_model()
% Returns a structure summarizing the contents of STACK, XXX, and META. 
% This function is often used to build an analysis summary cache file so 
% that the performance and coefficients of many models may be compared
% quickly without having to open each file each time.
%
% Iterates through the model, pulls out parameters that were fit, and 
% stores them in a structure for comparison of different fits.

global STACK XXX META;

if length(STACK) > length(XXX)
    recalc_xxx(1);
end

s = [];

% Extract summary of XXX
s.cellid = XXX{1}.cellid;
s.training_set = XXX{1}.training_set;
s.test_set = XXX{1}.test_set;
s.score_train_corr = XXX{end}.score_train_corr;
s.score_test_corr = XXX{end}.score_test_corr;
s.score_train_mse = XXX{end}.score_train_mse;
s.score_test_mse = XXX{end}.score_test_mse;

% Extract summary of STACK
s.n_free_params = length(pack_fittables(STACK));    
s.fits = [];
for ii = 1:length(STACK)
    if ~isfield(STACK{ii}, 'fit_fields')
        continue;
    end
    ff = STACK{ii}.fit_fields;
    for jj = 1:length(ff)
        s.fits.(STACK{ii}.name).(ff{jj}) = STACK{ii}.(ff{jj});
    end
end

% Extract summary of META
s.modelname   = META.modelname;
s.modelfile   = META.modelfile;
s.modelpath   = META.modelpath;
s.fit_time    = META.fit_time;
s.fitter      = META.fitter;
s.exit_code   = META.exit_code;
