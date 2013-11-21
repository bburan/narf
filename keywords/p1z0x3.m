function p1z0x3 ()

global MODULES STACK XXX;

% Cache the current fit parameters
savefitparms=cell(length(STACK),1);
for ii=1:length(STACK),
    if isfield(STACK{ii}{1},'fit_fields'),
        savefitparms{ii}=STACK{ii}{1}.fit_fields;
        STACK{ii}{1}.fit_fields={};
    end
end

append_module(MODULES.pole_zeros.mdl(...
                struct('fit_fields', {{'poles', 'gains', 'delays', 'y_offset'}}, ...
                       'n_poles', 1, ...
                       'n_zeros', 0)));

fit04a(); pop_module(); % Fit 1st only
    
% Now cache first P1Z0 params
savefitparms{end+1}=STACK{end}{1}.fit_fields;
STACK{end}{1}.fit_fields = {};

% Rename output
STACK{end}{1}.output = 'stim1';
calc_xxx(length(STACK)); 

% Add a second p1z0
append_module(MODULES.pole_zeros.mdl(...
                struct('fit_fields', {{'poles','gains', 'delays'}}, ...
                       'output', 'stim2', ...
                       'n_poles', 1, ...
                       'n_zeros', 0)));

append_module(MODULES.sum_fields.mdl(struct('inputs', {{'stim1', 'stim2'}})));

fit04a(); pop_module();  % Fit 2nd only
pop_module(); 

% Cache 2nd P1Z0 params
savefitparms{end+1}=STACK{end}{1}.fit_fields;
STACK{end}{1}.fit_fields = {};

% Add a third p1z0
append_module(MODULES.pole_zeros.mdl(...
                struct('fit_fields', {{'poles','gains', 'delays'}}, ...
                       'output', 'stim3', ...
                       'n_poles', 1, ...
                       'n_zeros', 0)));

append_module(MODULES.sum_fields.mdl(struct('inputs', {{'stim1', 'stim2', 'stim3'}})));

fit04a(); pop_module(); % Fit 3rd only

for ii=1:length(STACK)-2
    if isfield(STACK{ii}{1},'fit_fields'),
        STACK{ii}{1}.fit_fields=savefitparms{ii};
    end
end
      
fit04a();

end