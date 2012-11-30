function strf=tor_cellxc_demo(cellid)

%
% find the first TORC file for this cell in sCellFiles
%
dbopen(1);
sql=['SELECT * FROM sCellFile WHERE cellid="',cellid,'" AND runclassid=1'];

cellfiledata=mysql(sql);

ii=1;
rasterfs=100;

[parms,perf]=dbReadData(cellfiledata(ii).rawid);

parmfile=[cellfiledata(ii).stimpath cellfiledata(ii).stimfile];
spikefile=[cellfiledata(ii).path cellfiledata(ii).respfile];

%
% load the response PSTH
%
options=[];
options.channel=cellfiledata(ii).channum;
options.unit=cellfiledata(ii).unit;
options.rasterfs=rasterfs;
options.includeprestim=1;
options.tag_masks={'Reference'};

[resp,tags]=loadspikeraster(spikefile,options);

% convert to PSTH
resp=nanmean(resp,2);
resp=resp(:);

%
% load the stimulus
%
[stim,stimparam]=loadstimfrombaphy(parmfile,[],[],'specgram',rasterfs,24,0,1);
stim=stim(:,:)';

%
% choose fit algorithm and set various parameters
%
params=[];
params.altcore='cdcore';     % boosting, aka coordinate descent
%params.altcore='xccorefet';  % Theunissen et al 2001 algorithm

% min and max time lags in bins.  Min is typically zero
params.maxlag=[0 12];

% some other hyperparameters
params.resampcount=12;
params.sfscount=10;
params.sfsstep=3;


%
% calculate the STRF
%
strf=cellxcdataloaded(stim,resp,params);


%
% plot the STRF
%
figure

subplot(2,1,1);
plotastrf(strf(1).h,0,stimparam.ff,rasterfs);
title([cellid ' - STRF, no interpolation']);

subplot(2,1,2);
plotastrf(strf(1).h,2,stimparam.ff,rasterfs);
title([cellid ' - STRF, interpolated']);

