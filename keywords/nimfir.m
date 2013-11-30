function nimfir ()

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
append_module(MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
                                'baseline', 0,...
                                'fit_fields', {{'coefs','baseline'}})));

fit05(); pop_module();
                            
% Save the params for later, EXCEPT FOR THE Y_OFFSET!
savefitparms{end+1}=STACK{end}{1}.fit_fields;
STACK{end}{1}.fit_fields = {};
calc_xxx(length(STACK));

% Add normalization and a simple zsoft, then fit
append_module(MODULES.normalize_channels);
append_module(MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
                                              'phi', [1 1 0 0], ...
                                              'nlfn', @nl_softzero)));
fit05(); pop_module();

% Change the outputs to be named stim1
STACK{end-2}{1}.output = 'stim1'; % FIR
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

append_module(MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
                                'baseline', 0,...
                                'input', 'stim', 'output', 'stim2',...
                                'fit_fields', {{'coefs','baseline'}})));

append_module(MODULES.normalize_channels.mdl(struct('input', 'stim2', ...
                                                    'output', 'stim2')));
                                                
append_module(MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
                                              'phi', [-1 1 0 0], ...
                                              'input_stim', 'stim2', ...
                                              'output', 'stim2', ...
                                              'nlfn', @nl_softzero)));
                                              
append_module(MODULES.sum_fields.mdl(struct('inputs', {{'stim1', 'stim2'}})));

fit05(); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now fit everything together

for ii=1:6 % WONKY WONKY MAGIC NUMBER
    if isfield(STACK{ii}{1},'fit_fields'),
        STACK{ii}{1}.fit_fields=savefitparms{ii};
    end
end

fit05();

end