% function voltdepct2()
% 2nd-order volterra expansion followed by depct2 (not to be
% confused with depct2volt, which applies these in the opposite
% order)
%
function voltdepct2()

global MODULES XXX STACK;
global ANESTHETIZED_MOUSE

append_module(MODULES.add_nth_order_terms);

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

if ANESTHETIZED_MOUSE==1,
    disp('dep2: Using longer depression time constant IC for mouse data.');
    dep_tau_norm=100;
    fir_num_coefs=20;
    stop_exp=3.0;
else
    dep_tau_norm=100;
    fir_num_coefs=12;
    stop_exp=2.0;
end

append_module(MODULES.depression_filter_bank.mdl(...
    struct('strength', [1.5 10], ...
           'tau',      [10 20], ...
           'tau_norm',dep_tau_norm,...
           'crosstalk',20,...
           'num_channels', 1, ...
           'fit_fields', {{'strength', 'tau', 'crosstalk'}})));

meanresp = nanmean(flatten_field(XXX{end}.dat,XXX{end}.training_set,'respavg'));
append_module(MODULES.fir_filter.mdl(struct('num_coefs', fir_num_coefs, ...
                            'baseline',meanresp,...
                            'fit_fields', {{'coefs','baseline'}})));

fitSubstack([],10^-(stop_exp));
fitSubstack(length(STACK)-2,10^-(stop_exp+0.5));
