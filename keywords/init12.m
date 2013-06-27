function init12 ()
% Initializes all FIR filters via reverse correlation.

global STACK XXX;

[~, mod_idxs] = find_modules(STACK, 'fir_filter');

% Initialize coefs to all ones
for ii = 1:length(mod_idxs) % For each match
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
            h=zeros(size(strf(1).h));
            h(ii,:)=strf(1).h(ii,:);
            if ii==2,
                h(ii,:)=-h(ii,:);
            end
            STACK{idx}{jj}.coefs = h;
            STACK{idx}{jj}.num_dims = size(h, 1);
            
        end
    end
end

calc_xxx(mod_idxs{1});