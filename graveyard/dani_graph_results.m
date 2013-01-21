% Load all models in savepath, rebuild them, extract correlation scores
global NARF_PATH NARF_SAVED_MODELS_PATH XXX STACK;
savepath = NARF_SAVED_MODELS_PATH;
D = dir([savepath filesep '*.mat']);
sss = sprintf('\n\nCELLID    \tTRAIN CORR\tTRAINED ON\n');
for ii = 1:length(D),
    f = [savepath filesep D(ii).name];
    
    load_model_stack(f);
    recalc_xxx(1);
    
    if length(XXX{1}.training_set) > 1
        ts = 'all';
    else
        ts = XXX{1}.training_set{1};
    end
    
    c = STACK{end-1}.phi;
    
    sss = [sss sprintf('%s\t%f\t%s [%f %f %f %f]\n', ...
        XXX{1}.cellid, XXX{end}.score_train_corr, ts, ...
        c(1), c(2), c(3), c(4))];
end

fprintf('%s', sss);

open_narf_gui();
