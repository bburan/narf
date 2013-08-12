function update_metrics(batch, cellids, modelnames)

global STACK XXX META MODULES;

n_models = length(modelnames);
n_cellids = length(cellids);
stacks = cell(n_models, n_cellids);
metas = cell(n_models, n_cellids);
x0s = cell(n_models, n_cellids);

dbopen;
for mi = 1:n_models
    model = modelnames{mi};
    for ci = 1:n_cellids
        cellid = cellids{ci};
        
        sql = ['SELECT * FROM NarfResults WHERE batch=' num2str(batch) ''];
        sql = [sql ' AND cellid="' cellid '"'];
        sql = [sql ' AND modelname="' model '"'];
        results = mysql(sql);  
                
        fprintf('Updating Model [%d/%d]\n', (mi-1)*n_cellids+ci, n_models*n_cellids);        
        if isempty(results)
            stacks{mi, ci} = {};
            metas{mi, ci} = {};
            x0s{mi, ci} = {};
            fprintf('Skipping because it doesn''t exist.\n');
            continue;
        elseif (length(results) > 1)
            error('Multiple DB Results found: for %s/%s\n', cellid, model);
        end
        
        load_model(char(results(1).modelpath));
        
        % -------------------
        
%         append_module(MODULES.neural_statistics);        
% 
%         calc_xxx(1);
%         
%         fns = XXX{1}.training_set;
%         avgs = [];
%         for ii = 1:length(fns)
%             sf = fns{ii};
%             avgs(ii) = XXX{end}.dat.(sf).min_dist;
%         end        
%         META.metric_est_selfdist = sum(avgs); % Mathematically, this isn't true (distances do not always add), but I can't think of a fair way to compare stimuli from different data sets!
%         
%         fns = XXX{1}.test_set;
%         avgs = [];
%         for ii = 1:length(fns)
%             sf = fns{ii};
%             avgs(ii) = XXX{end}.dat.(sf).min_dist;
%         end        
%         META.metric_val_selfdist = sum(avgs); % Mathematically, this isn't true (distances do not always add), but I can't think of a fair way to compare stimuli from different data sets!
%         
        % ---------------------
        calc_xxx(1);   
        append_module(MODULES.load_stim_resps_from_baphy.mdl(...
                           struct('raw_resp_fs', 10000, ...
                                  'raw_stim_fs', 10, ...
                                  'include_prestim', true, ...
                                  'stimulus_format', 'envelope', ...
                                  'output_stim', 'stim10', ...
                                  'output_stim_time', 'stim10time', ...
                                  'output_resp', 'resp10000', ...
                                  'output_resp_time', 'resp10000time', ...
                                  'output_respavg', 'respavg'))); 
        append_module(MODULES.inter_spike_intervals.mdl(...
                           struct('input', 'resp10000', ...
                                  'time', 'resp10000time', ...
                                  'output', 'resp_ISIs')));

        append_module(MODULES.bayesian_likelihood);   % Default: Exponential
        STACK{end}{1}.probdist = 'exponential';
        STACK{end}{1}.probcutoff = 0.001; % Refractory period assumption
        calc_xxx(length(XXX)-3);
        META.metric_est_nlogl_exp = XXX{end}.score_train_nlogl; 
        META.metric_est_bic_exp =  XXX{end}.score_train_bic;
        META.metric_est_autocorr_exp =  XXX{end}.score_train_autocorr;
        META.metric_val_nlogl_exp = XXX{end}.score_test_nlogl; 
        META.metric_val_bic_exp =  XXX{end}.score_test_bic;
        META.metric_val_autocorr_exp =  XXX{end}.score_test_autocorr;
               
        STACK{end}{1}.probdist = 'lognormal';
        STACK{end}{1}.probcutoff = 0;
        calc_xxx(length(XXX)-1);        
        META.metric_est_nlogl_lognormal = XXX{end}.score_train_nlogl; 
        META.metric_est_bic_lognormal =  XXX{end}.score_train_bic;
        META.metric_est_autocorr_lognormal =  XXX{end}.score_train_autocorr;
        META.metric_val_nlogl_lognormal = XXX{end}.score_test_nlogl; 
        META.metric_val_bic_lognormal =  XXX{end}.score_test_bic;
        META.metric_val_autocorr_lognormal =  XXX{end}.score_test_autocorr;

        STACK{end}{1}.probdist = 'gamma';
        STACK{end}{1}.probcutoff = 0;
        calc_xxx(length(XXX)-1);                       
        META.metric_est_nlogl_gamma = XXX{end}.score_train_nlogl; 
        META.metric_est_bic_gamma =  XXX{end}.score_train_bic;
        META.metric_est_autocorr_gamma =  XXX{end}.score_train_autocorr;
        META.metric_val_nlogl_gamma = XXX{end}.score_test_nlogl; 
        META.metric_val_bic_gamma =  XXX{end}.score_test_bic;
        META.metric_val_autocorr_gamma =  XXX{end}.score_test_autocorr;

        STACK{end}{1}.probdist = 'inversegaussian';
        STACK{end}{1}.probcutoff = 0;
        calc_xxx(length(XXX)-1);                       
        META.metric_est_nlogl_inversegaussian = XXX{end}.score_train_nlogl; 
        META.metric_est_bic_inversegaussian =  XXX{end}.score_train_bic;
        META.metric_est_autocorr_inversegaussian =  XXX{end}.score_train_autocorr;
        META.metric_val_nlogl_inversegaussian = XXX{end}.score_test_nlogl; 
        META.metric_val_bic_inversegaussian =  XXX{end}.score_test_bic;
        META.metric_val_autocorr_inversegaussian =  XXX{end}.score_test_autocorr;
        
        STACK{end}{1}.probdist = 'weibull';
        STACK{end}{1}.probcutoff = 0;
        calc_xxx(length(XXX)-1);                       
        META.metric_est_nlogl_weibull = XXX{end}.score_train_nlogl; 
        META.metric_est_bic_weibull =  XXX{end}.score_train_bic;
        META.metric_est_autocorr_weibull =  XXX{end}.score_train_autocorr;
        META.metric_val_nlogl_weibull = XXX{end}.score_test_nlogl; 
        META.metric_val_bic_weibull =  XXX{end}.score_test_bic;
        META.metric_val_autocorr_weibull =  XXX{end}.score_test_autocorr;
                
        % ---------------------   
        
        mods = find_modules(STACK, 'fir_filter');
        mods = [mods{:}]; %% Flatten out parameter sets
        vals = cellfun(@(x) sparsity_metric(x.coefs), mods);
        META.metric_sparsity = sum(vals);
        vals = cellfun(@(x) smoothness_metric(x.coefs), mods);
        META.metric_smoothness = sum(vals);
        
        % ---------------------
        
        save_model(META.modelpath, STACK, XXX, META);
      
    end
end