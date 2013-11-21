function p1z0x2 ()

% A fairly intelligent keyword that:
%    1. Masks out all fit_fields except this module
%    2. Looks up the stack until a normalizer is found (usually after the
%       compressor nonlinearity and uses THAT as the STIM input. 
%    3. If this is the first pole_zero module, it fits the y_offset. 
%    4. Fits all parameters of p1. 
%    5. Restores all fit fields.
%    6. Fits all modules again. 
 
% (The idea is then you could chain them together and try BLAH_BLAH_BLAH and get better and better fits as you add more dimensions)
% Should need only 3 parameters: delay, gain, and pole (decay rate)

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