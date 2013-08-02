function ifn0()

global MODULES;
global STACK XXX;

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

append_module(MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
                                     'fit_fields', {{'coefs','baseline'}},...
                                     'output','stim1')));

append_module(MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
                                     'fit_fields', {{'coefs','baseline'}},...
                                     'output','stim2')));

[~, mod_idxs] = find_modules(STACK, 'fir_filter');

ii=1; % only initialize first fir filter.  Others kept at zero.

    for jj = 1:length(STACK{mod_idxs{ii}}) % For each paramset
        if isfield(STACK{mod_idxs{ii}}{jj}, 'fit_fields') && ...
                any(strcmp('coefs', STACK{mod_idxs{ii}}{jj}.fit_fields))
            mm = STACK{mod_idxs{ii}}{jj};
            idx = mod_idxs{ii};
            stim = [];
            resp = [];
            for kk = 1:length(XXX{idx}.training_set),
                f = XXX{idx}.training_set{kk};
                tstim=XXX{idx}.dat.(f).(mm.input);
                [Tx,Sx,Cx] = size(tstim);
                tstim=reshape(tstim, Tx*Sx, Cx);
                stim = cat(1, stim, tstim);
                resp = cat(1, resp, XXX{idx}.dat.(f).(mm.init_fit_sig)(:));
            end
            %resp = resp(:);
            %[Tx,Sx] = size(resp);
            %stim=reshape(permute(stim, [1 3 2]), Tx, numel(stim) / Tx);
            params = [];
            params.altcore     = 'cdcore';  % Either 'cdcore' or 'xccorefet'
            params.maxlag      = mm.num_coefs - 1;
            params.resampcount = 10;
            params.sfscount    = 5;
            params.sfsstep     = 3;
            strf = cellxcdataloaded(stim, resp, params);
            h=strf(1).h;
            STACK{idx}{jj}.coefs = h;
            STACK{idx}{jj}.num_dims = size(h, 1);
            
        end
    end

calc_xxx(mod_idxs{1});

% default inputs are already stim1, stim2
% default output already stim
append_module(MODULES.int_fire_neuron.mdl(struct(...
    'fit_fields', {{'Vrest','V0','gL'}},...
    'rectify_inputs',1)));
