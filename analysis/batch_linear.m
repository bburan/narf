% A batch file which tries to fit a linear model to multiple cellids
% Works only with SPN stim files

narf_set_path;
global NARF_PATH STACK XXX;

savepath = [NARF_PATH filesep 'saved_models'];

cellids = {'por026b-b2', ...
    'por026b-a1', ...
    'por025a-b1', ...
    'por024a-a1', ...
    'por024a-b1', ...
    'por023b-a1', ...
    'por022a-a1', ...
    'por022b-a1', ...
    'por020b-c1', ...
    'por019b-a1', ...
    'por019c-a1', ...
    'por018b-c1', ...
    'por017a-a1', ...
    'por016d-a1', ...
    'por016e-a1'};

% Function names to execute.
models = {'linear_fit_spn', ...
          'linear_fit_spn_stephen', ...
          'linear_fit_spn_exp', ...
          'linear_fit_spn_poly', ...
          'linear_fit_spn_heaviside', ...
          'linear_fit_spn_sigmoid', ...
          'linear_fit_spn_volterra', ...
          'linear_fit_spn_stephen_poly', ...
          'linear_fit_spn_inhib_excit', ...
          };

% UNCOMMENT FOR QUICK TESTING OF SCRIPT
models = {'linear_fit_spn_stephen'};
cellids = {'por025a-b1'};
      
scores = zeros(length(cellids), length(models) * 2);  

score_x_labels = cellids;
score_y_labels = cell(1, length(models) * 2);

for ci = 1:length(cellids),
    cid = cellids{ci};
    
    % Pick out the best training and test sets
    [cfd, cids, cellfileids] = dbgetscellfile('cellid', cid, ...
                                              'runclass','SPN');
    len = length(cfd);
    train_set = {};
    train_set_spikes = 0;
    test_set = {};
    test_set_reps = 0;
    
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
             cfd(i).spikes > train_set_spikes)
            train_set_spikes = cfd(i).spikes;
            train_set{1} = cfd(i).stimfile;
        end
    end
    
    % Fit using each of the models
    for mi = 1:length(models),
        m = models{mi};
    
        % Build the model and train it 
        fn = str2func(m);
        fn(cid, train_set);
        
        % Add the test set AFTER the training has occured, so that it is
        % not needlessly computed by default by the system. Recompute so
        % that the training and test scores are updated.
        XXX{1}.test_set = test_set;
        recalc_xxx(1); 
        
        % Save the file name 
        filename = sprintf('%s/%s__%f__%f__%s.mat', ...
            savepath, cid, ...
            XXX{end}.score_train_corr, ...
            XXX{end}.score_test_corr, ...
            m);
        
        save_model_stack(filename, STACK, XXX);
        
        % Save the test/train scores
        scores(ci, 2*mi-1) = XXX{end}.score_train_corr;
        scores(ci, 2*mi)  = XXX{end}.score_test_corr;
        
        score_y_labels{2*mi-1} = m;
        score_y_labels{2*mi} = sprintf('%s_test', m);
    end
end

% Plot the results
figure;
imagesc(scores');
%set(gca,'YDir','normal');
set(gca,'XTick', 1:length(score_x_labels));
set(gca,'YTick', 1:length(score_y_labels));
set(gca,'XTickLabel', score_x_labels);
set(gca,'YTickLabel', score_y_labels);
title('Scores');
axis tight;

% Build a list of the files in the 'modules/' directory
mdls = scan_directory_for_modules([NARF_PATH filesep 'modules/']);

pf = figure('Menubar','figure', 'Resize','off', ...
            'Units','pixels', 'Position', [20 50 1300 1100]);
narf_modelpane(pf, mdls); 