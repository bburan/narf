% A batch file which tries to fit a linear model to multiple cellids

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

for cid = cellids(:)', cid = cid{1};
    % Pick out the best training and test sets
    [cfd, cellids, cellfileids] = dbgetscellfile('cellid', cid);
    
    len = length(cfd);
    training_set = {};
    test_set = {};
    test_set_reps = 0;
    % parms = cell(1, len);
    % perfs = cell(1, len);
    
    % Load parms, perfs, if you need them
    % [parms{i}, perfs{i}] = dbReadData(cfd(i).rawid);

    % Select set with the most repetitions as test set.
    for i = 1:len;
        if (isfield(cfd(i), 'repcount') && cfd(i).repcount > test_set_reps)
           test_set_reps = dbget('gDataRaw', cfd(i).rawid, 'reps');
           test_set{1} = cfd(i).stimfile;
        end
    end

    % Train on every other passive trial by default
    for i = 1:len;
        if (~isequal(cfd(i).stimfile, test_set{1}) & ...
             isequal(cfd(i).behavior, 'passive'))
            training_set{end+1} = cfd(i).stimfile;
        end
    end

    train = {};
    test = {};
    
    % Fit using Stephen's boosting
    linear_fit_spn_stephen(cid, train);
    
    % TODO: Test the model on a different data set
    filename = sprintf('linear_fit_spn_stephen__%s__%f', cid, XXX{end}.score_corr);
    save_model_stack(filename, STACK)
    
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

