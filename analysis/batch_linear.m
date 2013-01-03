% A batch file which tries to fit a linear model to multiple cellids
% Works only with SPN stim files

narf_set_path;
global NARF_PATH STACK XXX;

cellids = {'por026b-b2', ...
    'por026b-a1', ...
    'por025a-b1', ...
    'por024b-c2', ...
    'por024a-a1', ...
    'por024a-b1', ...
    'por023b-a1', ...
    'por022a-a1', ...
    'por022b-a1', ...
    'por021a-b1', ...
    'por021b-b1', ...
    'por020b-c1', ...
    'por019b-a1', ...
    'por019c-a1', ...
    'por018b-c1', ...
    'por017a-a1', ...
    'por017c-1' , ...
    'por016d-a1', ...
    'por016e-a1'};

for ii = 1:length(cellids),
    cid = cellids(ii);
    
    % Pick out the best training and test sets
    [cfd, cellids, cellfileids] = dbgetscellfile('cellid', cid, ...
                                                 'runclass','SPN');
    len = length(cfd);
    training_set = {};
    training_set_spikes = 0;
    test_set = {};
    test_set_reps = 0;
    
    % Load parms, perfs, if you need them
    % parms = cell(1, len);
    % perfs = cell(1, len);
    % [parms{ii}, perfs{ii}] = dbReadData(cfd(ii).rawid);
    
    % Select the SPN sample with the most repetitions as test set.
    for i = 1:len;
        if (isfield(cfd(i), 'repcount') && cfd(i).repcount > test_set_reps)
           test_set_reps = cfd(i).repcount; 
           test_set{1} = cfd(i).stimfile;
        end
    end

    % Select the SPN with the most spikes as as the training set
    % (It also cannot be in the test set)
    for i = 1:len;
        if (~isequal(cfd(i).stimfile, test_set{1}) & ...
             cfd(i).spikes > training_set_spikes)
            training_set_spikes = cfd(i).spikes;
            training_set{1} = cfd(i).stimfile;
        end
    end
    
    % Fit using Stephen's boosting
    linear_fit_spn_stephen(cid, train);
    filename = sprintf('linear_fit_spn_stephen__%s__%f.mat', ...
                       cid, XXX{end}.score_corr);
    save_model_stack(filename, STACK);
    
%     % Fit using Matlab's lsqcurvefit
%     linear_fit_spn(cid);
%     filename = sprintf('linear_fit_spn__%s__%f', cid, XXX{end}.score_corr);
%     save_model_stack(filename, STACK)
%     
%     % Fit using Matlab's lsqcurvefit
%     linear_fit_spn_exp(cid);
%     filename = sprintf('linear_fit_spn_exp__%s__%f', cid, XXX{end}.score_corr);
%     save_model_stack(filename, STACK)
%     
%     % Fit using polynomial
%     linear_fit_spn_poly(cid);
%     filename = sprintf('linear_fit_spn_poly__%s__%f', cid, XXX{end}.score_corr);
%     save_model_stack(filename, STACK)
    
%     % Fit using volterra
    
end

% Build a list of the files in the 'modules/' directory
%mdls = scan_directory_for_modules([NARF_PATH filesep 'modules/']);

