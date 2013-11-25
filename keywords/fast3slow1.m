function fast3slow1 ()

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
append_module(MODULES.pole_zeros.mdl(...
                struct('fit_fields', {{'poles', 'zeros', 'gains', 'delays', 'y_offset'}}, ...
                       'n_poles', 3, ...
                       'n_zeros', 2)));                                      

% Set it to be manually be about the right filter 
                   
fit04a(); pop_module();
    
% Now cache those new params, EXCEPT FOR THE Y_OFFSET!
savefitparms{end+1}={{'poles', 'zeros', 'gains', 'delays'}};
STACK{end}{1}.fit_fields = {};
STACK{end}{1}.output = 'stim1';
calc_xxx(length(STACK)); 

% Add a second p1z0 (With a fittable y_offset)
append_module(MODULES.pole_zeros.mdl(...
                struct('fit_fields', {{'poles','gains', 'delays', 'y_offset'}}, ...
                       'output', 'stim2', ...
                       'n_poles', 1, ...
                       'n_zeros', 0)));
            
% Set it to be manually inhibitory and slow
STACK{end}{1}.poles = [-2];
STACK{end}{1}.gain = [-20];

append_module(MODULES.sum_fields.mdl(struct('inputs', {{'stim1', 'stim2'}})));

fit04a(); pop_module();

for ii=1:length(STACK)-2
    if isfield(STACK{ii}{1},'fit_fields'),
        STACK{ii}{1}.fit_fields=savefitparms{ii};
    end
end

% OPTIONAL TODO: Sum the two y_offsets together here before final fit.

fit04b();
end