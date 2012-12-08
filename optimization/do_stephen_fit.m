function coefs = do_stephen_fit(mdl, x)
% Uses Stephen's model fitting code.
%     mod_idx    the depth in the stack of stephen's FIR filter 
% Assumes all data has been prepared.

%global STACK XXX;
% mdl = STACK{mod_idx};
% x = xxx{mod_idx};

% Flatten each stimulus and response into a vector
mystim = {};
myresp = {};
fns = fieldnames(x.dat);
for ii = 1:length(fns)
    sf = fns{ii};
    
    mystim{ii} = permute(x.dat.(sf).ds_stim, [2 1 3]);
    dims = size(mystim{ii});
    mystim{ii} = reshape(mystim{ii}, dims(1)*dims(2), dims(3)); 
    
    myresp{ii} = x.dat.(sf).raw_respavg;
    myresp{ii} = myresp{ii}';
    myresp{ii} = myresp{ii}(:);    
end

% Concatenate all those together. The number of experimental dimensions
% must be equal for a concatenation to work; you can't mix a single-channel
% TORC and a multi-channel SPN. 
stim = [];
resp = [];
for ii = length(mystim)
    stim = cat(2, stim, mystim{ii}); % TODO: use a better algorithm and preallocate
    resp = cat(2, resp, myresp{ii});
end

% choose fit algorithm and set various parameters
params = [];
params.altcore     = mdl.altcore;  % Either 'cdcore' or 'xccorefet'
params.maxlag      = [0 mdl.maxlag];
params.resampcount = mdl.resampcount;
params.sfscount    = mdl.sfscount;
params.sfsstep     = mdl.sfsstep;

% Make an STRF
strf = cellxcdataloaded(stim,resp,params);

% Return the coefficients
coefs = strf(1).h;

end