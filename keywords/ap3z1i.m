function ap3z1i ()

global MODULES XXX STACK;

% Init the coefs to have the right dimensionality for 'stim'
x = XXX{end};
fns = fieldnames(x.dat);
sf = fns{1};
[T, S, num_chans] = size(x.dat.(sf).stim);                  

fir15;   %adds normalization and fir filter
fitSubstack([],10^-3);
meangain=sum(STACK{end}{1}.coefs,2);
meangain(meangain>=0)=0.6;
meangain(meangain<0)=-0.6;
pop_module;  % remove fir filter module, leave normalization

append_module(MODULES.pole_zeros.mdl(...
                struct('fit_fields', {{'gains', 'poles', 'zeros', 'delays'}}, ...
                       'n_poles', 3, ...
                       'n_zeros', 1, ...                       
                       'n_inputs', num_chans)));
STACK{end}{1}.gains(:)=meangain(:);
STACK{end}{1}.delays(:)=10;

append_module(MODULES.normalize_channels);


end