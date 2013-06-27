function isi_test()
% Just for development. Trying to figure out how bayesian log-likelihood
% can be used as an objective function

baphy_set_path;
narf_set_path; 

global MODULES;
MODULES = scan_directory_for_modules();

cellid = 'por025a-b1';
training_set = {'por025a07_p_SPN'};
test_set = {'por025a06_p_SPN'};

% Define how groups of modules should be combined to make many models
STACK = {}; 
mm{1}.isi100 = {MODULES.load_stim_resps_from_baphy.mdl(...
                                     struct('raw_resp_fs', 10000, ...
                                            'raw_stim_fs', 100,...
                                            'stimulus_format', 'envelope')), 
                MODULES.inter_spike_intervals};  
                                        
mm{2} = module_groups('log2');
mm{3} = module_groups('firb');
mm{4} = module_groups('nonl');
mm{5} = module_groups('fminlsq');
mm{6} = module_groups('mse');

[~, modelnames] = module_combinations(mm);

fit_models(mm, cellid, training_set, test_set);
               
open_narf_gui();