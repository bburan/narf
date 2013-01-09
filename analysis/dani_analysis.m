% Fit a model to five stimulus files to initialize model params.
% Fit a specific model parameter again using each file individually.
% Save the fitted params and plot the results.

narf_set_path;
global NARF_PATH STACK XXX;
mdls = scan_directory_for_modules([NARF_PATH filesep 'modules']);

savepath = [NARF_PATH filesep 'saved_models'];

files_to_analyze = [...
    {'por022a-a1', 'por022a07_p_SPN', 	'por022a09_a_TSP', 	'por022a12_p_SPN', 	'por022a13_a_TSP',	'por022a14_p_SPN'};
    {'por022a-c1', 'por022a07_p_SPN'	'por022a09_a_TSP', 	'por022a12_p_SPN',	'por022a13_a_TSP',	'por022a14_p_SPN'};
    %{'por022b-a1', 'por022b08_p_SPN', 	'por022b10_a_TSP', 	'por022b11_p_SPN', 	'por022b12_a_TSP', 	'por022b13_p_SPN'}; % por022b12_a_TSP has weird resp size!
    %{'por022b-a2', 'por022b08_p_SPN',	'por022b10_a_TSP', 	'por022b11_p_SPN', 	'por022b12_a_TSP', 	'por022b13_p_SPN'}; 
    %{'por023a-a1', 'por023a06_p_SPN', 	'por023a07_a_TSP', 	'por023a08_p_SPN', 	'por023a09_a_TSP',	'por023a10_p_SPN'}; % por023a08_p_SPN is not yet sorted!
    %{'por023a-b1', 'por023a06_p_SPN', 	'por023a07_a_TSP', 	'por023a08_p_SPN', 	'por023a09_a_TSP',	'por023a10_p_SPN'};
    %{'por023a-c1', 'por023a06_p_SPN', 	'por023a07_a_TSP', 	'por023a08_p_SPN', 	'por023a09_a_TSP',	'por023a10_p_SPN'};
    {'por023b-a1', 'por023b12_p_SPN', 	'por023b13_a_TSP', 	'por023b15_p_SPN',	'por023b16_a_TSP', 	'por023b18_p_SPN'};
    {'por023b-b1', 'por023b12_p_SPN', 	'por023b13_a_TSP', 	'por023b15_p_SPN',	'por023b16_a_TSP', 	'por023b18_p_SPN'};
    {'por023b-d1', 'por023b12_p_SPN', 	'por023b13_a_TSP', 	'por023b15_p_SPN',	'por023b16_a_TSP', 	'por023b18_p_SPN'};
    {'por023b-d2', 'por023b12_p_SPN', 	'por023b13_a_TSP', 	'por023b15_p_SPN',	'por023b16_a_TSP', 	'por023b18_p_SPN'};
    %{'por024b-b1', 'por024b03_p_SPN', 	'por024b04_a_TSP', 	'por024b05_p_SPN', 	'por024b07_p_TSP', 	'por024b08_p_SPN'}; % por024b07_p_TSP has weird size
    %{'por024b-c1', 'por024b03_p_SPN', 	'por024b04_a_TSP', 	'por024b05_p_SPN', 	'por024b07_p_TSP', 	'por024b08_p_SPN'};
    {'por025a-b1', 'por025a09_p_SPN', 	'por025a10_a_TSP', 	'por025a11_p_SPN', 	'por025a12_a_TSP', 	'por025a13_p_SPN'};
    {'por025a-c1', 'por025a09_p_SPN', 	'por025a10_a_TSP', 	'por025a11_p_SPN', 	'por025a12_a_TSP', 	'por025a13_p_SPN'};
    {'por025a-c2', 'por025a09_p_SPN', 	'por025a10_a_TSP', 	'por025a11_p_SPN', 	'por025a12_a_TSP', 	'por025a13_p_SPN'};
    {'por025a-d1', 'por025a09_p_SPN', 	'por025a10_a_TSP', 	'por025a11_p_SPN', 	'por025a12_a_TSP', 	'por025a13_p_SPN'};
    %{'por025c-b1', 'por025c04_p_SPN',	'por025c05_a_TSP', 	'por025c06_p_SPN', 	'por025c07_a_TSP', 	'por025c08_a_SPN'}; % por025c07_a_TSP has weird size
    %{'por025c-c1', 'por025c04_p_SPN',	'por025c05_a_TSP', 	'por025c06_p_SPN', 	'por025c07_a_TSP', 	'por025c08_a_SPN'};
    %{'por025c-c2', 'por025c04_p_SPN',	'por025c05_a_TSP', 	'por025c06_p_SPN', 	'por025c07_a_TSP', 	'por025c08_a_SPN'};
    {'por026a-b1', 'por026a08_p_SPN', 	'por026a09_a_TSP', 	'por026a10_p_SPN', 	'por026a11_a_TSP',	'por026a12_p_SPN'};
    {'por026a-d1', 'por026a08_p_SPN', 	'por026a09_a_TSP', 	'por026a10_p_SPN', 	'por026a11_a_TSP',	'por026a12_p_SPN'}];

% -------------------------------------------------------------------------
% Define the model
raster_fs = 200;
filter_length = 10;
n_channels = 2;

XXX = {};
XXX{1} = [];

STACK = {};
STACK{1} = mdls.load_stim_resps_from_baphy.mdl(...
                struct('raw_resp_fs', raster_fs, ...
                       'raw_stim_fs', raster_fs,...
                       'stimulus_format','envelope', ...
                       'stimulus_channel_count', n_channels));
STACK{2} = mdls.normalize_channels;
STACK{3} = mdls.fir_filter.mdl(struct('num_dims', n_channels, ...
                                      'num_coefs', filter_length));
STACK{4} = mdls.correlation;
% -------------------------------------------------------------------------

[M, N] = size(files_to_analyze);  

% For each cell file
for ci = 1:M,  
    cellid = files_to_analyze{ci, 1};
    respfiles = {files_to_analyze{ci, 2:6}};
    % Build the model, train it on everything, and save 
    XXX{1}.cellid = cellid;  
    XXX{1}.training_set = respfiles;
    XXX{1}.test_set = {};
    STACK{3}.fit_fields = {'coefs'};
    fit_with_lsqcurvefit();
    filename = sprintf('%s/%s_all.mat', savepath, cellid);
    save_model_stack(filename, STACK, XXX);
    
    % Now train the model again on each respfile individually, with
    % possibly different free model parameters
    %STACK{3}.fit_fields = {};
    for rfi = 2:6
        rf = files_to_analyze{ci, rfi};
        XXX = {};
        XXX{1} = [];
        XXX{1}.training_set = {rf};
        XXX{1}.cellid = cellid;  
        XXX{1}.test_set = {};
        fit_with_lsqcurvefit();
        filename = sprintf('%s/%s_%s.mat', savepath, cellid, rf);
        save_model_stack(filename, STACK, XXX);
    end
end

% Open up a display so you can view the results
% pf = figure('Menubar','figure', 'Resize','off', ...
%             'Units','pixels', 'Position', [20 50 1300 min(170*length(STACK)+40, 1100)]);
% narf_modelpane(pf, mdls); 

% Now load each module again, rebuild and extract their correlation scores
D = dir([savepath filesep '*.mat']);
sss = sprintf('\n\nCELLID\tTRAIN CORR\tTRAINED ON\n');
for ii = 1:length(D),
    f = [savepath filesep D(ii).name];
    
    load_model_stack(f);
    recalc_xxx(1);
    
    if length(XXX{1}.training_set) > 1
        ts = 'all';
    else
        ts = XXX{1}.training_set{1};
    end
    
    sss = [sss sprintf('%s\t%f\t%s\n', XXX{1}.cellid, XXX{end}.score_train_corr, ts)];
end

fprintf('%s', sss);




