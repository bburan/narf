% cellxcdemo.m
%
% example usage of modelgenresponse.m and cellxcnodb.m
%
% path requirements:
%   ... /cellxc/
%   ... /mutils/
%   ... /gen/
%   ... /db/
%
% created SVD 3/11/03
%

% make sure parameters are being fully reset
clear params xcparms

% generate model data
params.stimfile='/auto/k5/david/data/modelhi/test.natrev.andros.30_pix.imsm';
params.outfile='/tmp/natrev.andros.30_pix.resp.mat';
params.cellid='model';
params.threshold=30;
params.ori=0;
params.sf=3;
params.modelnoise=1;
params.simple=0;

modelgenresponse(params);


% estimate the model kernel
xcparms.stimfiles={params.stimfile};
xcparms.respfiles={params.outfile};
xcparms.resploadcmd='respload';
xcparms.resploadparms={'',0,1,1};
xcparms.stimloadcmd='loadimfile';
xcparms.stimloadparms={32,0,16};

%xcparms.kernfmt='space';
xcparms.kernfmt='pfft';
%xcparms.kernfmt='pfft+4';
if strcmp(xcparms.kernfmt,'space'),
   xcparms.stimfiltercmd='';
   xcparms.stimfilterparms={};
elseif strcmp(xcparms.kernfmt,'pfft'),
   xcparms.stimfiltercmd='movpower';
   xcparms.stimfilterparms={0,0,1,1,50};
else
   xcparms.stimfiltercmd='movphasesep';
   xcparms.stimfilterparms={0,0,1,1,0};
end

% this may help but runs slower
xcparms.altfit='xcfit2';

xcparms.outfile='/tmp/modelout.mat';
xcparms.showres=1;
xcparms.sfsstep=7;
xcparms.sfscount=10;
xcparms.sffiltsigma=5;
xcparms.smoothtime=0;
xcparms.fitfrac=0;
xcparms.predfrac=0.1;

cellxcnodb(xcparms);

return




%
% real data natrev demo
%

% fill out xcparms from above. then add additional fields:
xcparms.cellid='r0210a';
% runclassid specifies type of data (0=review, 2=natrev, etc)
xcparms.runclassid=2;
xcparms.resploadparms={'',1,1,1};
xcparms.outfile=['/tmp/',xcparms.cellid,'.mat'];

[cellfiles, times, xcparms] = cellfiletimes(xcparms.cellid, xcparms)

cellxcnodb(xcparms);

return




%
%  video demo
%

% estimate the model kernel
xcparms.stimfiles={'/auto/k5/hayden/van_data/2003-11-23/v0098.video4.007.imsm'};
xcparms.respfiles={'/auto/k5/hayden/van_data/2003-11-23/v0098.video4.007.resp.mat'};
xcparms.resploadcmd='respload';
xcparms.resploadparms={'',0,1,1};
xcparms.stimloadcmd='loadimfile';
xcparms.stimloadparms={32,0,16};

xcparms.stimfiltercmd='movpower';
xcparms.stimfilterparms={0,0,1,1,0};
xcparms.kernfmt='pfft';

% this may help but runs slower
xcparms.altfit='xcfit2';

xcparms.outfile='/tmp/modelout.mat';
xcparms.showres=1;
xcparms.sfsstep=8;
xcparms.sfscount=20;
xcparms.sffiltsigma=5;
xcparms.smoothtime=0;
xcparms.fitfrac=0;
xcparms.predfrac=0.1;

cellxcnodb(xcparms);

