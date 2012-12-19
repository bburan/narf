% Fits an STRF from the preprocessed stimulus, not the envelope.
% Uses an elliptical filter

narf_set_path;

global NARF_PATH STACK XXX;

% Build a list of the files in the 'modules/' directory
mdls = scan_directory_for_modules([NARF_PATH filesep 'modules/']);

% XXX defines the initial data
XXX = {};
XXX{1}.cellid = 'por024a-a1';
XXX{1}.training_set = {'por024a19_p_SPN'};
XXX{1}.test_set = {};

% Define the model structure, using defaults for unspecified fields.
n_channels = 2;
filter_length = 20;
raster_fs = 200;

STACK = {};
STACK{1} = mdls.load_stim_resps_from_baphy;

% Autodetect the low and high frequencies 
recalc_xxx(1);
sf = XXX{1}.training_set{1}; %% (Assumes 1st training set is SPN)
hit_idx = arrayfun(@(x)strcmp(x.stimfile, sf), XXX{2}.cfd);
[parm,perf]=dbReadData(XXX{2}.cfd(hit_idx).rawid);
low_freqs = parm.Ref_LowFreq;
high_freqs = parm.Ref_HighFreq;
    
STACK{2} = mdls.elliptic_bandpass_filter_bank.mdl(struct('num_channels', n_channels, ...
                                                 'low_freqs', low_freqs, ...
                                                 'high_freqs', high_freqs));                                 
STACK{3} = mdls.downsample_with_fn;
STACK{4} = mdls.fir_filter.mdl(struct('num_dims', n_channels, ...
                                      'num_coefs', filter_length));

STACK{5} = mdls.correlation.mdl(struct('input1', 'stim', ...
                                       'input2', 'respavg'));

% NOTE: Must use double braces {{ }} because of struct()'s default behavior
% of unfortunately stripping away the first layer of cell references.

% Compute the entire stack once
recalc_xxx(1);

m = [];
m.num_dims = n_channels;
m.altcore = 'cdcore';
m.maxlag = filter_length - 1;
m.resampcount = filter_length - 1;  
m.sfscount = 10;
m.sfsstep = 3;
m.rasterfs = raster_fs;
m.coefs = zeros(n_channels, filter_length);
m.stimfield = 'stim'; % Names of the fields to use as stimulus and response
m.respfield = 'resp';

% Run Stephen's fitting routine 
STACK{4}.coefs = do_stephen_fit(m, XXX{2});
ss = size(STACK{4}.coefs);
STACK{4}.num_coefs = ss(1);
STACK{4}.num_dims = ss(2);

% Recompute now just the FIR filte
recalc_xxx(4); 

% Finally, display the GUI for easy tweaking and viewing of the best result
pf = figure('Menubar','figure', 'Resize','off', ...
            'Units','pixels', 'Position', [20 50 1300 900]);
narf_modelpane(pf, mdls); 
