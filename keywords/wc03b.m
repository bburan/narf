function wc03b()
% Weight all input channels, producing 1 output channel
% Works on 'stim' by default. 
global MODULES STACK XXX;
signal = 'stim';
fns = fieldnames(XXX{end}.dat);
dschans=8;
n_output_chans = 3;

% Fit an FIR filter first and use its principal components
n_input_chans = size(XXX{end}.dat.(fns{1}).(signal), 3);
mapbins=ceil((1:n_input_chans)./n_input_chans.*dschans);
Vds=zeros(n_input_chans,dschans);
for ii=1:dschans,
    Vds(mapbins==ii,ii)=1./sum(mapbins==ii);
end

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));
append_module(MODULES.weight_channels.mdl(struct('weights', Vds)));

firtemp2();

coefs = STACK{end}{1}.coefs;
coefs=Vds*coefs;
[U,S,V]=svd(coefs);

pop_module();
pop_module();

append_module(MODULES.weight_channels.mdl(...
       struct('weights', U(:, 1:n_output_chans), ...
              'fit_fields', {{'weights'}})));

