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

if 0,
   % generate model data
   params.cellid='model';
   params.freq=[3500 4500];
   params.fbw=[0.4 0.6];
   params.lat=[25 35];
   params.tdecay=[13 10];
   params.amp=[1 -1];
   params.tbincount=13;
   params.tfs=100;
   params.threshold=-100;
   params.stimloadcmd='loadspeech';
   params.stimloadparms={'specgram',100,32,0};
   params.trialcount=10;
   params.gaussnoisestd=3;
   
   respfiles={};
   
   params.stimfile='/afs/glue.umd.edu/department/isr/labs/nsl/projects/daqsc/data16/Coral/c033a05_p_sp1';
   params.outfile='/auto/k5/david/data/model/model.sp1.mat';
   modelA1(params);
   stimfiles{1}=params.stimfile;
   respfiles{1}=params.outfile;
   
   params.stimfile='/afs/glue.umd.edu/department/isr/labs/nsl/projects/daqsc/data16/Coral/c033a05_p_sp2';
   params.outfile='/auto/k5/david/data/model/model.sp2.mat';
   modelA1(params);
   stimfiles{2}=params.stimfile;
   respfiles{2}=params.outfile;
else
   stimfiles={'/afs/glue.umd.edu/department/isr/labs/nsl/projects/daqsc/data16/Coral/c033a05_p_sp1',...
              '/afs/glue.umd.edu/department/isr/labs/nsl/projects/daqsc/data16/Coral/c033a05_p_sp2'};
   respfiles={'/auto/k5/david/data/model/model.sp1.mat',...
              '/auto/k5/david/data/model/model.sp2.mat'};
   params.stimloadcmd='loadspeech';
   params.stimloadparms={'specgram',100,32,0};
end

disp('ONLY USING ONE FILE!!!');
stimfiles={stimfiles{1}};
respfiles={respfiles{1}};

% set up parameters for cellxcnodb
xcparms.respfiles=respfiles;
xcparms.stimfiles=stimfiles;
xcparms.resploadcmd='respload';
xcparms.resploadparms={'',0,1,1};

xcparms.stimloadcmd=params.stimloadcmd;
xcparms.stimloadparms=params.stimloadparms;
xcparms.stimfiltercmd='';
xcparms.stimfilterparms={};

xcparms.kernfmt='spect';


% cross val fit
xcparms.altfit='xcfit2';

if 0,
   xcparms.altcore='xccore';
   xcparms.minlag=-6;
   xcparms.sfscount=15;
else
   xcparms.altcore='cdcore';
   xcparms.minlag=0;
   xcparms.sfscount=30;
   %xcparms.stimfiltercmd='gsmooth';
   %xcparms.stimfilterparms={[0.5 1],0,0};
end

xcparms.outfile='/auto/k5/david/data/model/model.sp2.out.mat';
xcparms.showres=1;
xcparms.sfsstep=8;
xcparms.sffiltsigma=5;
xcparms.smoothtime=0;
xcparms.fitfrac=0;
xcparms.predfrac=0.1;
xcparms.decorrspace=4; % need full space-time decorr
xcparms.maxlag=12;
xcparms.resampcount=10;

% estimate the model kernel
cellxcnodb(xcparms);


figure(2);
clf

subplot(2,1,2);
zz=zload(xcparms.outfile);
strfest=zz.strf(2).h;
mm=max(abs(strfest(:)));
imagesc(strfest,[-mm mm]);
colormap('default');
axis xy;


load(xcparms.respfiles{1});
[stim,ff]=feval(params.stimloadcmd,params.stimfile,params.respstart,...
                params.respstop,params.stimloadparms{:});
psthtrue=kernpredict(strf,stim,1,0,1);
psthest=kernpredict(strfest,stim,1,0,1);

subplot(2,1,1);
mm=max(abs(strf(:)));
imagesc(strf,[-mm mm]);
colormap('default');
axis xy;

return
load(xcparms.respfiles{1});
[stim,ff]=feval(params.stimloadcmd,params.stimfile,1,0,params.stimloadparms{:});
psthtrue=kernpredict(strf,stim,1,0,1);
psthest=kernpredict(strfest,stim,1,0,1);
