clear all 
close all

% Fit a model to five stimulus files to initialize model params.
% Fit a specific model parameter again using each file individually.
% Save the fitted params and plot the results.

narf_set_path;
global NARF_PATH STACK XXX;
mdls = scan_directory_for_modules([NARF_PATH filesep 'modules']);


files_to_analyze = [...
    % OLD SINGLE-SPEAKER RECORDINGS
    %     {'por020a-a1',	 'por020a04_p_SPN',	'por020a08_a_TSP',	'por020a09_p_SPN',	'por020a10_a_TSP', 	'por020a11_p_SPN'};
    %    {'por020a-c1',  	 'por020a04_p_SPN',	'por020a08_a_TSP',  'por020a09_p_SPN',	'por020a10_a_TSP', 	'por020a11_p_SPN'};
    %     {'por020a-d1', 	 'por020a04_p_SPN',	'por020a08_a_TSP',  'por020a09_p_SPN',	'por020a10_a_TSP', 	'por020a11_p_SPN'};
    %     {'por019b-a1', 	 'por019b02_p_SPN',	'por019b05_a_TSP',  'por019b07_p_SPN', 	'por019b08_a_TSP', 	'por019b09_p_SPN'};
    %     {'por019b-a3', 	 'por019b02_p_SPN',	'por019b05_a_TSP',  'por019b07_p_SPN', 	'por019b08_a_TSP', 	'por019b09_p_SPN'};
    %     {'por019b-b1', 	 'por019b02_p_SPN',	'por019b05_a_TSP',  'por019b07_p_SPN',	'por019b08_a_TSP', 	'por019b09_p_SPN'};
    %     {'por019c-a1', 	 'por019c04_p_SPN', 'por019c06_a_TSP',  'por019c07_p_SPN', 	'por019c08_a_TSP', 	'por019c09_p_SPN'};
    %     {'por018b-a1',     'por018b04_p_SPN',	'por018b06_a_TSP',  'por018b07_p_SPN',	'por018b08_a_TSP',	'por018b09_p_SPN'};
    %     {'por018b-a2',     'por018b04_p_SPN',	'por018b06_a_TSP',  'por018b07_p_SPN',	'por018b08_a_TSP',	'por018b09_p_SPN'};
    %     {'por018b-c1',     'por018b04_p_SPN',	'por018b06_a_TSP',  'por018b07_p_SPN',	'por018b08_a_TSP',	'por018b09_p_SPN'};
    %     {'por018b-d1',     'por018b04_p_SPN',	'por018b06_a_TSP',  'por018b07_p_SPN',	'por018b08_a_TSP',	'por018b09_p_SPN'};
    %     {'por018c-a1',     'por018c05_p_SPN',	'por018c06_a_TSP',  'por018c08_p_SPN',	'por018c09_p_TSP',	'por018c10_p_SPN'};
    %     {'por018c-a2',     'por018c05_p_SPN',	'por018c06_a_TSP',  'por018c08_p_SPN',	'por018c09_p_TSP',	'por018c10_p_SPN'};
    %     {'por018c-c1',     'por018c05_p_SPN',	'por018c06_a_TSP',  'por018c08_p_SPN',	'por018c09_p_TSP',	'por018c10_p_SPN'};
    %     {'por018c-d1',  	 'por018c05_p_SPN',	'por018c06_a_TSP',  'por018c08_p_SPN',	'por018c09_p_TSP',	'por018c10_p_SPN'};
    
    % TWO-SPEAKER RECORDINGS
    
    %{'por022a-a1', 'por022a08_p_SPN', 	'por022a09_a_TSP', 	'por022a12_p_SPN', 	'por022a13_a_TSP',	'por022a14_p_SPN'};  
    %{'por022a-c1', 'por022a08_p_SPN'	'por022a09_a_TSP', 	'por022a12_p_SPN',	'por022a13_a_TSP',	'por022a14_p_SPN'};
    %{'por022b-a1', 'por022b08_p_SPN', 	'por022b10_a_TSP', 	'por022b11_p_SPN', 	'por022b12_a_TSP', 	'por022b13_p_SPN'}; % crash; por022b12_a_TSP has weird resp size! 
    %{'por022b-a2', 'por022b08_p_SPN',	'por022b10_a_TSP', 	'por022b11_p_SPN', 	'por022b12_a_TSP', 	'por022b13_p_SPN'}; % crash
    %{'por023a-a1', 'por023a06_p_SPN', 	'por023a07_a_TSP',  'por023a09_p_TSP',	'por023a10_p_SPN'};                     % por023a08_p_SPN is bad, skip; crummy STRF for this channel
    %{'por023a-b1', 'por023a06_p_SPN', 	'por023a07_a_TSP', 	'por023a09_p_TSP',	'por023a10_p_SPN'};                     
    %{'por023a-c1', 'por023a06_p_SPN', 	'por023a07_a_TSP', 	'por023a09_p_TSP',	'por023a10_p_SPN'};                     % really BAD fitting!
    %{'por023b-a1', 'por023b12_p_SPN', 	'por023b13_a_TSP', 	'por023b15_p_SPN',	'por023b16_a_TSP', 	'por023b18_p_SPN'};
    %{'por023b-b1', 'por023b12_p_SPN', 	'por023b13_a_TSP', 	'por023b15_p_SPN',	'por023b16_a_TSP', 	'por023b18_p_SPN'};
    %{'por023b-d1', 'por023b12_p_SPN', 	'por023b13_a_TSP', 	'por023b15_p_SPN',	'por023b16_a_TSP', 	'por023b18_p_SPN'}; 
    %{'por023b-d2', 'por023b12_p_SPN', 	'por023b13_a_TSP', 	'por023b15_p_SPN',	'por023b16_a_TSP', 	'por023b18_p_SPN'}; 
    %{'por024b-b1', 'por024b03_p_SPN', 	'por024b04_a_TSP', 	'por024b05_p_SPN', 	'por024b07_p_TSP', 	'por024b08_p_SPN'}; % crash; por024b07_p_TSP has weird size!
    %{'por024b-c1', 'por024b03_p_SPN', 	'por024b04_a_TSP', 	'por024b05_p_SPN', 	'por024b07_p_TSP', 	'por024b08_p_SPN'}; % crash
    %{'por025a-b1', 'por025a09_p_SPN', 	'por025a10_a_TSP', 	'por025a11_p_SPN', 	'por025a12_a_TSP', 	'por025a13_p_SPN'};
    %{'por025a-c1', 'por025a09_p_SPN', 	'por025a10_a_TSP', 	'por025a11_p_SPN', 	'por025a12_a_TSP', 	'por025a13_p_SPN'}; % 
    %{'por025a-c2', 'por025a09_p_SPN', 	'por025a10_a_TSP', 	'por025a11_p_SPN', 	'por025a12_a_TSP', 	'por025a13_p_SPN'}; % really BAD fitting!
    %{'por025a-d1', 'por025a09_p_SPN', 	'por025a10_a_TSP', 	'por025a11_p_SPN', 	'por025a12_a_TSP', 	'por025a13_p_SPN'};
    %{'por025c-b1', 'por025c04_p_SPN',	'por025c05_a_TSP', 	'por025c06_p_SPN', 	'por025c07_a_TSP', 	'por025c08_a_SPN'}; % crash; por025c07_a_TSP has weird size
    %{'por025c-c1', 'por025c04_p_SPN',	'por025c05_a_TSP', 	'por025c06_p_SPN', 	'por025c07_a_TSP', 	'por025c08_a_SPN'}; % crash
    %{'por025c-c2', 'por025c04_p_SPN',	'por025c05_a_TSP', 	'por025c06_p_SPN', 	'por025c07_a_TSP', 	'por025c08_a_SPN'}; % crash
    %{'por026a-b1', 'por026a08_p_SPN', 	'por026a09_a_TSP', 	'por026a10_p_SPN', 	'por026a11_a_TSP',	'por026a12_p_SPN'};
    %{'por026a-d1', 'por026a08_p_SPN', 	'por026a09_a_TSP', 	'por026a10_p_SPN', 	'por026a11_a_TSP',	'por026a12_p_SPN'};
    %{'por027a-a1', 'por027a06_p_SPN',	'por027a07_a_TSP',	'por027a08_p_SPN',	'por027a09_a_TSP',	'por027a13_p_SPN'}; % crash; Unit a1 lost in B2
    %{'por027a-b1', 'por027a06_p_SPN',	'por027a07_a_TSP',	'por027a08_p_SPN',	'por027a09_a_TSP',	'por027a13_p_SPN'}; % crash; Unit b1 lost in B2
    %{'por027a-b2', 'por027a06_p_SPN',	'por027a07_a_TSP',	'por027a08_p_SPN',	'por027a09_a_TSP',	'por027a13_p_SPN'}; % crash; Unit b2 lost in B2
    %{'por027a-c1', 'por027a06_p_SPN',	'por027a07_a_TSP',	'por027a08_p_SPN',	'por027a09_a_TSP',	'por027a13_p_SPN'};
    {'por027b-b1', 'por027b06_p_SPN', 	'por027b08_a_TSP',	'por027b11_p_SPN',	'por027b14_a_TSP',	'por027b16_p_SPN'}; % crash
    %{'por027b-b2', 'por027b06_p_SPN', 	'por027b08_a_TSP',	'por027b11_p_SPN',	'por027b14_a_TSP',	'por027b16_p_SPN'}; % crash
    %{'por027b-c1', 'por027b06_p_SPN', 	'por027b08_a_TSP',	'por027b11_p_SPN',	'por027b14_a_TSP',	'por027b16_p_SPN'}; % crash
    %{'por028b-a1', 'por028b02_p_SPN',	'por028b05_a_TSP',	'por028b06_p_SPN',	'por028b10_a_TSP',	'por028b12_p_SPN'}; % crash
    %{'por028b-b1', 'por028b02_p_SPN',	'por028b05_a_TSP',	'por028b06_p_SPN',	'por028b10_a_TSP',	'por028b12_p_SPN'}; % crash
    %{'por028b-b2', 'por028b02_p_SPN',	'por028b05_a_TSP',	'por028b06_p_SPN',	'por028b10_a_TSP',	'por028b12_p_SPN'}; % crash
    %{'por028b-c1', 'por028b02_p_SPN',	'por028b05_a_TSP',	'por028b06_p_SPN',	'por028b10_a_TSP',	'por028b12_p_SPN'}; % crash
    %{'por028b-d1', 'por028b02_p_SPN',	'por028b05_a_TSP',	'por028b06_p_SPN',	'por028b10_a_TSP',	'por028b12_p_SPN'}; % crash
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
                       'include_prestim',0,...
                       'stimulus_format','envelope', ...
                       'stimulus_channel_count', n_channels));
%STACK{2} = mdls.depression_filter_bank             
STACK{2} = mdls.normalize_channels; % Uncomment to not consider depression
STACK{3} = mdls.fir_filter.mdl(struct('num_dims', n_channels, ...
                                      'num_coefs', filter_length));
%STACK{5} = mdls.nonlinearity.mdl(struct('phi', [0.05 0.5 0.05 0], ...
%                                        'nlfn', @sigmoidal));
STACK{4} = mdls.correlation;

% -------------------------------------------------------------------------

[M, N] = size(files_to_analyze);

% For each cell file
for ci = 1:M,  
    cellid = files_to_analyze{ci, 1};
    savepath = [NARF_PATH filesep 'saved_models' filesep cellid];
    if ~exist(savepath,'dir'),
        mkdir(fileparts(savepath),cellid);
    end
    filecount=length({files_to_analyze{ci,:}});
    respfiles = {files_to_analyze{ci, 2:filecount}};
    
    % Re-initialize the data to evaluate
    XXX{1}.cellid = cellid;  
    XXX{1}.training_set = respfiles;
    XXX{1}.test_set = {};
    
     % Initial fit with stephen's second method
    recalc_xxx(1);
    stim = [];
    resp = [];
    for ii = 1:length(XXX{1}.training_set),
        sf = XXX{1}.training_set{ii};
        stim = cat(1, stim, XXX{3}.dat.(sf).stim);
        % Make them wider than they need to be
        [M N P] = size(XXX{3}.dat.(sf).resp);
        temp = nan * zeros(M,N,4); % 4 is a magic number. It must be larger than any expected resp files dimension 3
        temp(~isnan(XXX{3}.dat.(sf).resp)) = XXX{3}.dat.(sf).resp(~isnan(XXX{3}.dat.(sf).resp));
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
    STACK{3}.coefs = strf(1).h;
    recalc_xxx(3);  % Recompute now from the FIR filter onward
    
    filename = sprintf('%s/%s_all_xccorefet.mat', savepath, cellid);
    save_model_stack(filename, STACK, XXX);
    
    % re-fit fir filter using least-square
    STACK{3}.fit_fields = {'coefs'};
    fit_with_lsqcurvefit();
    
    % add output NL to stack
    linpredmin=min(XXX{4}.dat.(sf).stim(:));
    linpredmax=max(XXX{4}.dat.(sf).stim(:));
    phi1=mean(std(XXX{4}.dat.(sf).stim(:)))-std(XXX{4}.dat.(sf).stim(:))*2;
    phi2=0.001;
    phi3=max(resp(:))/6;
    STACK{4} = mdls.nonlinearity.mdl(struct('phi', [phi1 phi2 phi3 0], ...
                                            'nlfn', @sigmoidal));
    STACK{5} = mdls.correlation;
    
    STACK{3}.fit_fields = {};
    STACK{4}.fit_fields = {'phi'};
    fit_with_lsqcurvefit();
    phi_init=STACK{4}.phi;
    
    filename = sprintf('%s/%s_all.mat', savepath, cellid);
    save_model_stack(filename, STACK, XXX);
    
    %open_narf_gui();
    
    % Now train the model again on each respfile individually, holding the
    % depression and FIR coefficients free and letting the nonlinearity slide
    nlbins=1000;
    nlmatrix=zeros(nlbins,5);
    filecount=length({files_to_analyze{ci,:}});
    for rfi = 2:filecount,
        rf = files_to_analyze{ci, rfi};
        XXX = {};
        XXX{1} = [];
        XXX{1}.training_set = {rf};
        XXX{1}.cellid = cellid;  
        XXX{1}.test_set = {};
        STACK{4}.phi = phi_init;  % Reset phi each time to avoid minima
        fit_with_lsqcurvefit();
        filename = sprintf('%s/%s_%s.mat', savepath, cellid, rf);
        save_model_stack(filename, STACK, XXX);
        
        % save some results for plotting
        x=linspace(linpredmin,linpredmax,nlbins)';
        nlmatrix(:,rfi-1)=sigmoidal(STACK{4}.phi,x);
        
    end
end

filename = sprintf('%s/%s_all.mat', savepath, cellid);
load_model_stack(filename);
open_narf_gui();
set(gcf,'name',[cellid ' : ' basename(filename)]);

figure;
plot(x,nlmatrix);
filecount=length({files_to_analyze{ci,:}});
hl=legend({files_to_analyze{ci, 2:filecount}},'Location','Northwest');
set(hl,'interpreter','none');
set(gcf,'name',[cellid ' : gain per file']);

