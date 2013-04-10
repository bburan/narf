function initrc ()
% Initializes all FIR filters via reverse correlation.

global STACK, XXX;

% TODO: Search through STACK and find all FIR filters with fircoefs in fit_fields

% TODO: init each FIR filter using reverse correlation

% Initialize coefs automatically if it's it has fit_fields
if any(strcmp('coefs', mm.fit_fields))
    stim = [];
    resp = [];
    for ii = 1:length(xxx{end}.training_set),
        f = xxx{end}.training_set{ii};
        tstim=xxx{end}.dat.(f).(mm.input);
        [Tx,Sx,Cx] = size(tstim);
        tstim=reshape(tstim, Tx*Sx, Cx);
        stim = cat(1, stim, tstim);
        resp = cat(1, resp, xxx{end}.dat.(f).(mm.init_fit_sig)(:));
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
    mm.coefs = strf(1).h;
    mm.num_dims = size(mm.coefs, 1);
    
end

