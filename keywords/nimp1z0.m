function nimp1z0 ()

global MODULES STACK XXX;

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

% Cache the current fit parameters
savefitparms=cell(length(STACK),1);
for ii=1:length(STACK),
    if isfield(STACK{ii}{1},'fit_fields'),
        savefitparms{ii}=STACK{ii}{1}.fit_fields;
        STACK{ii}{1}.fit_fields={};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHWAY 1 (EXCITATORY)
% Start by fitting the fast part of the filter
append_module(MODULES.pole_zeros.mdl(...
                struct('fit_fields', {{'poles', 'gains', 'delays'}}, ...
                       'n_poles', 1, ...
                       'n_zeros', 0)));                               
fit04a(); pop_module();

% Save the params for later, EXCEPT FOR THE Y_OFFSET!
savefitparms{end+1}={'poles', 'gains', 'delays'}; % Remove yoffset
STACK{end}{1}.fit_fields = {}; % Stop fitting poles
calc_xxx(length(STACK));

% Add normalization and a simple zsoft, then fit
append_module(MODULES.normalize_channels);
append_module(MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
                                              'phi', [1 1 0 0], ...
                                              'nlfn', @nl_softzero)));
fit04a(); pop_module();

% Change the outputs to be named stim1
STACK{end-2}{1}.output = 'stim1'; % polezero
STACK{end-1}{1}.output = 'stim1'; % Normalizer
STACK{end-1}{1}.input = 'stim1';
STACK{end}{1}.output = 'stim1'; % Nonlinearity
STACK{end}{1}.input = 'stim1';
savefitparms{end+1} = {}; % For the damn normalizer
savefitparms{end+1} = STACK{end}{1}.fit_fields; % Push 
STACK{end}{1}.fit_fields = {}; % Don't fit nonlinearity now
calc_xxx(length(STACK)-2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHWAY 2 (FAST INHIBITORY)

append_module(MODULES.pole_zeros.mdl(...
                struct('fit_fields', {{'poles', 'gains', 'delays'}}, ...
                       'input', 'stim', ...
                       'output', 'stim2', ...
                       'n_poles', 1, ...
                       'n_zeros', 0)));    

append_module(MODULES.normalize_channels.mdl(struct('input', 'stim2', ...
                                                    'output', 'stim2')));
append_module(MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
                                              'phi', [-1 1 0 0], ...
                                              'input_stim', 'stim2', ...
                                              'output', 'stim2', ...
                                              'nlfn', @nl_softzero)));
                                              
append_module(MODULES.sum_fields.mdl(struct('inputs', {{'stim1', 'stim2'}})));

fit04a(); pop_module(); pop_module();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHWAY 3 (Slow Inhibitory Dynamics)

% Add a second p1z0 (With a fittable y_offset)
append_module(MODULES.pole_zeros.mdl(...
                struct('fit_fields', {{'poles','gains', 'delays', 'y_offset'}}, ...
                       'output', 'stim3', ...
                       'n_poles', 1, ...
                       'n_zeros', 0)));
            
% Set it to be manually inhibitory and slow
STACK{end}{1}.poles = [-2]; % Nice and slow
STACK{end}{1}.gain = [-20];

append_module(MODULES.sum_fields.mdl(struct('inputs', {{'stim1', 'stim2', 'stim3'}})));

for ii=1:6 % WONKY WONKY MAGIC NUMBER
    if isfield(STACK{ii}{1},'fit_fields'),
        STACK{ii}{1}.fit_fields=savefitparms{ii};
    end
end

% Refit everything hard
fit04b();
end