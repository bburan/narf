function nim()

global MODULES STACK XXX;

NMMPATH = '/home/ivar/matlab/nmm/';
addpath(genpath(NMMPATH));

% TODO: Assumes 200Hz for now. 
dt = 0.005;
input_data = [];
input_data.stim = [];
input_data.zvector = [];
input_data.spiketimes_ms = {};

% USES THE TRAINING DATA ONLY OF COURSE
sf = XXX{end}.training_set;
for ii = 1:length(sf),
    sf = sf{ii};
    % TODO: Presently, this only supports a single channel
    tmp = XXX{end}.dat.(sf).stim(:,:,1);
    input_data.stim = cat(1, input_data.stim, tmp(:));
end
input_data.zvector = ones(size(input_data.stim));

% To create the input_data.spiketimes_ms, temporarily create modules
% to get the RESP file at 10KHz sampling rate as spike-times again
append_module(MODULES.load_stim_resps_from_baphy.mdl(...
                           struct('raw_resp_fs', 10000, ...
                                  'raw_stim_fs', 200,...
                                  'include_prestim', true, ...
                                  'stimulus_format', 'envelope'))); 

append_module(MODULES.inter_spike_intervals); 

input_data.spiketimes_ms = {};
for ii = 1:length(sf),
    for rr = 1:size(XXX{end}.dat.(sf).stim, 3)
        input_data.spiketimes_ms{rr} = [];       
        for ss = 1:size(XXX{end}.dat.(sf).stim, 2)
            
            t_offset = (ss-1) * (size(XXX{end}.dat.(sf).stim, 1) * dt);
            tmp = XXX{end}.dat.(sf).resp_spiketimes(:,ss,rr);
            idxs = tmp~=0;
            input_data.spiketimes_ms{rr} = ...
                    cat(1, input_data.spiketimes_ms{rr}, 1000*(t_offset + tmp(idxs)));
        end
    end
end

pop_module(); % Remove last modules
pop_module(); 

%% Parameters for fit
SR = 200;  % sampling frequency
tent_basis_spacing = 1; % represent stimulus filters using tent-bases with this spacing (in up-sampled time units)
nLags = 30; % number of time lags for estimating stimulus filters
frac = 2; % amount of up-sampling of stim resolution

% Param structure for 1-channel (temporal only)
params_stim = NMMcreate_stim_params( nLags, 1/SR, frac, tent_basis_spacing );

% Param structure for 2-channel stimulus
%params_stim2 = NMMcreate_stim_params( [nLags 2], 1/SR, frac, tent_basis_spacing );

% Shift response by in order to "make room" for beginning of filter -- particularly important PSC term
sh = 0;

%% Parse data into correct format
% Function I designed to parse SVD data structure. Embedded in this is the all-important 
% "create_time_embedding" function that generates design (X) matrix with stim-params

params_stim.frac = frac;

[Xtot spks msd rspks stim] = A1Tdesign(input_data, params_stim, 1 );

% spks0 = data.spiketimes_ms{n}'/1000;
    
Robs_tot = histc(spks,(0:1/frac:4*length(stim))/SR);
Robs_tot = shift_mat_zpad(Robs_tot,sh);

% Define one example of nested 10-fold x-val (randomly selected time points)
NT = size(Xtot,1); Nxv = ceil(NT*0.1);
rorder = randperm(NT);
training_set = rorder(1:(NT-Nxv));
xv_set = rorder((NT-Nxv+1):end);
% Nested within training set
NestTrain = rorder(1:NT-2*Nxv);
NestTest = rorder((NT-2*Nxv+1):(NT-Nxv));

Xtrain = Xtot(training_set,:);    Rtrain = Robs_tot(training_set);
XNtrain = Xtot(NestTrain,:);      RNtrain = Robs_tot(NestTrain);
Xtest = Xtot(xv_set,:);           Rtest = Robs_tot(xv_set);
XNtest = Xtot(NestTest,:);        RNtest = Robs_tot(NestTest);

%% Begin fitting
% Initialize model
params_reg = NMMcreate_reg_params( 'lambda_d2T', 100, 'lambda_L1', 0 );
params_reg.boundary_conds(1) = 0; % make zero-temporal boundary conditions (for smoothness reg)
fit0 = NMMCinitialize_model( params_stim, 1, {'lin'}, params_reg, 1 );
silent = 0;

% Fit GLM (filters)
fit0 = NMMCfit_filters( fit0, Rtrain, Xtrain, [], [], silent );

% Make rectified and refit
fit0r = fit0;  fit0r.mods(1).NLtype = 'threshlin';
fit0r = NMMCfit_filters( fit0r, Rtrain, Xtrain, [], [], silent );
% This hurts fit -- stick with linear

%% Adjust filter regularization for kicks here (using nested)
L1s = [0 1 4 10 20 40 100]; % see if L1 helps at all
L2s = [1 2 5 10 50 100 500 1000]; % Check range

nestlist = zeros(length(L1s),length(L2s)); 
for mm = 1:length(L1s)
  for nn = 1:length(L2s)
    test = NMMadjust_regularization( fit0, 1, 'lambda_d2T', L2s(nn), 'lambda_L1', L1s(mm) );
    test = NMMCfit_filters( test, RNtrain, XNtrain, [], [], 1 );
    %% test nested cross-validation data
    nestlist(mm,nn) = NMMCmodel_eval( test, RNtest, XNtest ); 
    fprintf( '  L1 %2d | L2 %4d: %f\n', L1s(mm), L2s(nn), nestlist(mm,nn) )
  end
end
[a bestL2indx] = max(max(nestlist));
[a bestL1indx] = max(nestlist(:,bestL2indx));
fprintf( 'Best regularization: L1 = %d  L2 = %d\n', L1s(bestL1indx), L2s(bestL2indx) )

fit0 = NMMadjust_regularization( fit0, 1, 'lambda_d2T', L2s(bestL2indx), 'lambda_L1', L1s(bestL1indx) );
fit0 = NMMCfit_filters( fit0, Rtrain, Xtrain, [], [], 1 );
% Wow -- ends up being pretty temporally quick, and little bumpy elsewhere. see below for more comments

% Try fitting spike-history term
fitH = NIMinit_spkhist( fit0, 12, 1, 6, 1 );
fitH = NMMCfit_filters( fitH, Rtrain, Xtrain, [], [], silent );
% note that first bin is not suppressed (RF violation?) -- 25 spikes within 0.5 ms...

% Ignore spike history for now -- check upstream_LNs
fit1 = NMMinitialize_upstreamNLs( fit0, Xtot, 1, 40 );
fit1 = NMMCfit_upstreamNLs( fit1, Rtrain, Xtrain, [], [], silent );
% this finds largely linear with saturation
% lets leave as linear term when fitting next one

%% Add delayed inhibition term
fit2 = NMMCadd_NLinput( fit0, 'threshlin', -1, 1, shift_mat_zpad( fit0.mods(1).filtK, 3 ) );
fit2 = NMMCfit_filters( fit2, Rtrain, Xtrain, [], [], silent );
% Note that it inherited the regualization from fit0, but might be better to be different

% Just check on second adjustments
nestlist = zeros(length(L1s),length(L2s)); 
for mm = 1:length(L1s)
  for nn = 1:length(L2s)
    test = NMMCadd_NLinput( fit0, 'threshlin', -1, 1, shift_mat_zpad( fit0.mods(1).filtK, 3 ) );
    test = NMMadjust_regularization( test, 2, 'lambda_d2T', L2s(nn), 'lambda_L1', L1s(mm) );
    test = NMMCfit_filters( test, RNtrain, XNtrain, [], [], 1 );
    %% test nested cross-validation data
    nestlist(mm,nn) = NMMCmodel_eval( test, RNtest, XNtest ); 
    fprintf( '  L1 %2d | L2 %4d: %f\n', L1s(mm), L2s(nn), nestlist(mm,nn) )
  end
end
[a bestL2indx] = max(max(nestlist));
[a bestL1indx] = max(nestlist(:,bestL2indx));
fprintf( 'Best regularization: L1 = %d  L2 = %d\n', L1s(bestL1indx), L2s(bestL2indx) )

fit2 = NMMCadd_NLinput( fit0, 'threshlin', -1, 1, shift_mat_zpad( fit0.mods(1).filtK, 3 ) );
fit2 = NMMadjust_regularization( fit2, 2, 'lambda_d2T', L2s(bestL2indx), 'lambda_L1', L1s(bestL1indx) );
fit2 = NMMCfit_filters( fit2, Rtrain, Xtrain, [], [], 1 );

% Initialize upstream nonlinearities and see if helps (alternating optimizations)
% note that it seems to help to optimize inhibition before making exc nonlinear in this case
fit2u = NMMinitialize_upstreamNLs( fit2, Xtot, 2, 10 );
fit2u = NMMCfit_alt( fit2u, Rtrain, Xtrain, [], [], silent );
% Note that leaving exc as linear is helpful here -- but then fit both
fit2u = NMMinitialize_upstreamNLs( fit2u, Xtot, 1, 10 );
fit2u = NMMCfit_alt( fit2u, Rtrain, Xtrain, [], [], silent );
% I didn't explore the regularzation for the nonlinearity, but its usually much less sensitive
% Nor did I optimize number of lags -- looks like even less might be helpful so it doesn't require
% such strong L1 regularization. Still, this is what exploratory fits look like.
% Nested cross-validation also might want to use a little more data

%displayNMM(fit2u);

% INSERT the trained model into the wrapper module for saving. 
append_module(MODULES.nim_wrapper.mdl(struct('nim_model', fit2u)));
