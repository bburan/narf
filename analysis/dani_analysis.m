% Fit a model to five stimulus files to initialize model params.
% Fit a specific model parameter again using each file individually.
% Save the fitted params and plot the results.

narf_set_path;
global NARF_PATH STACK XXX;
mdls = scan_directory_for_modules([NARF_PATH filesep 'modules']);

savepath = [NARF_PATH filesep 'saved_models'];

files_to_analyze = [...
    {'por022a-a1', 'por022a08_p_SPN', 	'por022a09_a_TSP', 	'por022a12_p_SPN', 	'por022a13_a_TSP',	'por022a14_p_SPN'}; % PROBLEM: 
    %{'por022a-c1', 'por022a08_p_SPN'	'por022a09_a_TSP', 	'por022a12_p_SPN',	'por022a13_a_TSP',	'por022a14_p_SPN'};
    %{'por022b-a1', 'por022b08_p_SPN', 	'por022b10_a_TSP', 	'por022b11_p_SPN', 	'por022b12_a_TSP', 	'por022b13_p_SPN'}; % por022b12_a_TSP has weird resp size! 
    %{'por022b-a2', 'por022b08_p_SPN',	'por022b10_a_TSP', 	'por022b11_p_SPN', 	'por022b12_a_TSP', 	'por022b13_p_SPN'}; 
    %{'por023a-a1', 'por023a06_p_SPN', 	'por023a07_a_TSP', 	'por023a08_p_SPN', 	'por023a09_a_TSP',	'por023a10_p_SPN'}; % por023a08_p_SPN is not yet sorted!
    %{'por023a-b1', 'por023a06_p_SPN', 	'por023a07_a_TSP', 	'por023a08_p_SPN', 	'por023a09_a_TSP',	'por023a10_p_SPN'};
    %{'por023a-c1', 'por023a06_p_SPN', 	'por023a07_a_TSP', 	'por023a08_p_SPN', 	'por023a09_a_TSP',	'por023a10_p_SPN'};
    %{'por023b-a1', 'por023b12_p_SPN', 	'por023b13_a_TSP', 	'por023b15_p_SPN',	'por023b16_a_TSP', 	'por023b18_p_SPN'};
    %{'por023b-b1', 'por023b12_p_SPN', 	'por023b13_a_TSP', 	'por023b15_p_SPN',	'por023b16_a_TSP', 	'por023b18_p_SPN'};
    %{'por023b-d1', 'por023b12_p_SPN', 	'por023b13_a_TSP', 	'por023b15_p_SPN',	'por023b16_a_TSP', 	'por023b18_p_SPN'};
    %{'por023b-d2', 'por023b12_p_SPN', 	'por023b13_a_TSP', 	'por023b15_p_SPN',	'por023b16_a_TSP', 	'por023b18_p_SPN'};
    %{'por024b-b1', 'por024b03_p_SPN', 	'por024b04_a_TSP', 	'por024b05_p_SPN', 	'por024b07_p_TSP', 	'por024b08_p_SPN'}; % por024b07_p_TSP has weird size
    %{'por024b-c1', 'por024b03_p_SPN', 	'por024b04_a_TSP', 	'por024b05_p_SPN', 	'por024b07_p_TSP', 	'por024b08_p_SPN'};
    %{'por025a-b1', 'por025a09_p_SPN', 	'por025a10_a_TSP', 	'por025a11_p_SPN', 	'por025a12_a_TSP', 	'por025a13_p_SPN'};
    %{'por025a-c1', 'por025a09_p_SPN', 	'por025a10_a_TSP', 	'por025a11_p_SPN', 	'por025a12_a_TSP', 	'por025a13_p_SPN'};
    %{'por025a-c2', 'por025a09_p_SPN', 	'por025a10_a_TSP', 	'por025a11_p_SPN', 	'por025a12_a_TSP', 	'por025a13_p_SPN'};
    %{'por025a-d1', 'por025a09_p_SPN', 	'por025a10_a_TSP', 	'por025a11_p_SPN', 	'por025a12_a_TSP', 	'por025a13_p_SPN'};
    %{'por025c-b1', 'por025c04_p_SPN',	'por025c05_a_TSP', 	'por025c06_p_SPN', 	'por025c07_a_TSP', 	'por025c08_a_SPN'}; % por025c07_a_TSP has weird size
    %{'por025c-c1', 'por025c04_p_SPN',	'por025c05_a_TSP', 	'por025c06_p_SPN', 	'por025c07_a_TSP', 	'por025c08_a_SPN'};
    %{'por025c-c2', 'por025c04_p_SPN',	'por025c05_a_TSP', 	'por025c06_p_SPN', 	'por025c07_a_TSP', 	'por025c08_a_SPN'};
    %{'por026a-b1', 'por026a08_p_SPN', 	'por026a09_a_TSP', 	'por026a10_p_SPN', 	'por026a11_a_TSP',	'por026a12_p_SPN'};
    %{'por026a-d1', 'por026a08_p_SPN', 	'por026a09_a_TSP', 	'por026a10_p_SPN', 	'por026a11_a_TSP',	'por026a12_p_SPN'};
    ];

% -------------------------------------------------------------------------
% Define the model
raster_fs = 100;
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
STACK{2} = mdls.depression_filter_bank             
STACK{2} = mdls.normalize_channels; % Uncomment to not consider depression
STACK{3} = mdls.normalize_channels;
STACK{4} = mdls.fir_filter.mdl(struct('num_dims', n_channels, ...
                                      'num_coefs', filter_length));
STACK{5} = mdls.nonlinearity.mdl(struct('phi', [0.05 0.5 0.05 0], ...
                                        'nlfn', @sigmoidal));
STACK{6} = mdls.correlation;
% -------------------------------------------------------------------------

[M, N] = size(files_to_analyze);

% For each cell file
for ci = 1:M,  
    cellid = files_to_analyze{ci, 1};
    respfiles = {files_to_analyze{ci, 2:6}};
    
    % Re-initialize the data to evaluate
    XXX{1}.cellid = cellid;  
    XXX{1}.training_set = respfiles;
    XXX{1}.test_set = {};
    
    %STACK{2}.fit_fields = {'tau'};
    STACK{4}.fit_fields = {'coefs'};
    
    % Initial fit with stephen's second method
    recalc_xxx(1);
    stim = [];
    resp = [];
    for ii = 1:length(XXX{1}.training_set),
        sf = XXX{1}.training_set{ii};
        stim = cat(1, stim, XXX{4}.dat.(sf).stim);
        % Make them wider than they need to be
        [M N P] = size(XXX{4}.dat.(sf).resp);
        temp = nan * zeros(M,N,4); % 4 is a magic number. It must be larger than any expected resp files dimension 3
        temp(~isnan(XXX{4}.dat.(sf).resp)) = XXX{4}.dat.(sf).resp(~isnan(XXX{4}.dat.(sf).resp));
        resp = cat(1, resp, temp);
    end
    resp=nanmean(resp,3);
    resp=resp(:);
    [T,S] = size(resp);
    stim=reshape(permute(stim, [1 3 2]), T, numel(stim) / T);
    params = [];
    params.altcore     = 'xccorefet';  % Either 'cdcore' or 'xccorefet'
    params.maxlag      = filter_length - 1;
    params.resampcount = filter_length - 1;
    params.sfscount    = 10;
    params.sfsstep     = 3;
    strf = cellxcdataloaded(stim, resp, params);
    STACK{4}.coefs = strf(1).h;
    
    STACK{5}.fit_fields = {'phi'};
    
    filename = sprintf('%s/%s_all_xccorefet.mat', savepath, cellid);
    save_model_stack(filename, STACK, XXX);
    
    % Ironically, this seems to work worse than nothing at all:
    % Try to guess some better initial conditions for the optimization
    %     % FIR DEFAULT: Cols 2,3 start as ones, everything else zero
    %     STACK{4}.coefs = zeros(size(STACK{4}.coefs));
    %     STACK{4}.coefs(:, 2:3) = ones(size(STACK{4}.coefs(:, 2:3)));
    %     
    %     % For SIGMA, 
    %     recalc_xxx(1);
    %     rr = flatten_field(XXX{5}.dat, XXX{1}.training_set, 'respavg');
    %     ss = flatten_field(XXX{5}.dat, XXX{1}.training_set, 'stim');
    %     mu = mean(ss);
    %     sigma = std(ss);
    %     R = corrcoef(ss, rr);
    %     amp = 2 * mean(rr) * sign(R(2,1));
    %     offset = 0;
    %     STACK{5}.phi = [mu sigma amp offset];
         
    fit_with_lsqcurvefit();
    filename = sprintf('%s/%s_all.mat', savepath, cellid);
    save_model_stack(filename, STACK, XXX);
    
    % Now train the model again on each respfile individually, holding the
    % depression and FIR coefficients free and letting the nonlinearity slide
    STACK{2}.fit_fields = {};
    STACK{4}.fit_fields = {};
    STACK{5}.fit_fields = {'phi'};
    for rfi = 2:6
        rf = files_to_analyze{ci, rfi};
        XXX = {};
        XXX{1} = [];
        XXX{1}.training_set = {rf};
        XXX{1}.cellid = cellid;  
        XXX{1}.test_set = {};
        STACK{5}.phi = [0.05 0.5 0.05 0];  % Reset phi each time to avoid minima
        fit_with_lsqcurvefit();
        filename = sprintf('%s/%s_%s.mat', savepath, cellid, rf);
        save_model_stack(filename, STACK, XXX);
    end
end

open_narf_gui();


