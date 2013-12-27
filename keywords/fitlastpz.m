function fitlastpz(input_channels)

global STACK MODULES;

% The output is a normalized, weighted linear combination of the above. 
append_module(MODULES.concatenate_channels.mdl(struct('inputs', {input_channels})));
append_module(MODULES.normalize_channels);
wc01();
fitSubstack(length(STACK),10^-5); % Fit just the weights

fit04a(); pop_module(); % Pop NMSE
pop_module();  % Pop wc01
pop_module();  % Pop normalize
pop_module();  % Pop concat channels
