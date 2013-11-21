function p1z0x2 ()

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

fit04a(); pop_module();
    
% Now cache those new params
savefitparms{end+1}=STACK{end}{1}.fit_fields;
STACK{end}{1}.fit_fields = {};
STACK{end}{1}.output = 'stim1';
calc_xxx(length(STACK)); 

% Add a second p1z0
append_module(MODULES.pole_zeros.mdl(...
                struct('fit_fields', {{'poles','gains', 'delays'}}, ...
                       'output', 'stim2', ...
                       'n_poles', 1, ...
                       'n_zeros', 0)));

append_module(MODULES.sum_fields.mdl(struct('inputs', {{'stim1', 'stim2'}})));

fit04a(); pop_module();

for ii=1:length(STACK)-2
    if isfield(STACK{ii}{1},'fit_fields'),
        STACK{ii}{1}.fit_fields=savefitparms{ii};
    end
end
      
fit04a();
end