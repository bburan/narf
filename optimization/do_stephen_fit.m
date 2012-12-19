function coefs = do_stephen_fit(mdl, x)
% Uses Stephen's model fitting code to find filter coefficients

fns = fieldnames(x.dat);

% Concatenate all input files together. The number of experimental dimensions
% must be equal for a concatenation to work; you can't mix a single-channel
% TORC and a multi-channel SPN. Sorry. 
stim = [];
resp = [];

% Concatenate the stimulus of every stim file
% TODO: Complain to stephen that overlapping could be a problem
for sf = fieldnames(x.dat)', sf = sf{1};
%     ARBITRARY DIMENSION STUB ----------------------------
%     % Get the stimulus matrix and column selector/accessor
%     M = x.dat.(sf).(mdl.stimfield);
%     M_selcols = x.dat.(sf).([mdl.stimfield '_selcols']);
%         
%     % Each channel must become a very long single stimulus. 
%     for chan_idx = 1:mdl.num_dims
%         tmp = M(:, M_selcols('chan', chan_idx));
%         tmpstim(:, chan_idx) = tmp(:);
%     end
%     ARBITARY_DIMENSION_STUB------------------------------
    stim = cat(1, stim, x.dat.(sf).(mdl.stimfield));
    resp = cat(1, resp, x.dat.(sf).(mdl.respfield));
end

% choose fit algorithm and set various parameters
params = [];
params.altcore     = mdl.altcore;  % Either 'cdcore' or 'xccorefet'
params.maxlag      = [0 mdl.maxlag];
params.resampcount = mdl.resampcount;
params.sfscount    = mdl.sfscount;
params.sfsstep     = mdl.sfsstep;

% Adjust the data's format to match Stephen's interface
resp=nanmean(resp,3);
resp=resp(:);
[T,S] = size(resp);
stim=reshape(permute(stim, [1 3 2]), T, numel(stim) / T);

% Make an STRF
strf = cellxcdataloaded(stim, resp, params);

% Return the coefficients
coefs = strf(1).h;

end