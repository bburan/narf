function dep2ifn()

global MODULES STACK;

% intialize excitatory input same as standard dep2
dep2();

% code to split pos and negative components of IC to separate
% inputs.  currently not being used
allcoefs=STACK{end}{1}.coefs;
poscoefs=allcoefs.*(allcoefs>0);
negcoefs=-allcoefs.*(allcoefs<0);
%STACK{end}{1}.coefs=poscoefs;
%STACK{end}{1}.baseline=0;

% re-route output to 'stim1', which means excitatory synapse
STACK{end}{1}.output='stim1';
calc_xxx(2);

% initialize inhibitory synapse to be all zero.  better
% alternative? One option is to set it to the negative components
% of linear fit (STACK{end}{1}.coefs=negcoefs), but this seems to
% hurt
append_module(MODULES.fir_filter.mdl(struct('num_coefs', 12,...
                                     'fit_fields', {{'coefs','baseline'}},...
                                     'output','stim2')));

% default inputs to IFN are already stim1, stim2
% default output already stim
append_module(MODULES.int_fire_neuron.mdl(struct(...
    'fit_fields', {{'Vrest','V0','gL'}},...
    'rectify_inputs',1)));
fitSubstack();

