% function [hf,stephf,movformat,kern]=wgxconline(rawids,fmt);
%
% created SVD 5/13/02
%
% for WG: use fmt='parm'
% for NR: use fmt='sfft' for phase separated fourier domain
%         use fmt='wav' for wavelet domain
%
function [hf,stephf,movformat,kern]=wgxconline(rawids,fmt);

% rawfileids gives entries into gDataRaw for data sets to analyze
if ~exist('rawids','var'),
   rawfileids=[596,597,598];
   %rawfileids=[592];
end
if ~exist('fmt','var'),
   fmt='parm';
end

% cellfileids give entries into sCellFile for data sets to analyze
cellfileids=dbraw2scellfile(rawids,fmt);

dbopen;
BATQUEUEID=[];

for ii=1:2,
   figure(ii);
end
drawnow;

% set up all the parameters for cellxc.m:   
sfileid='(';
for ii=1:length(cellfileids),
   sfileid=[sfileid,num2str(cellfileids(ii)),','];
end
sfileid(end)=')';
sql=['SELECT * FROM sCellFile WHERE id in ',sfileid];
cellfiledata=mysql(sql);

stimfiles={};
respfiles={};
resplens=zeros(length(cellfiledata),1);
for ii=1:length(cellfiledata),
   stimfiles{ii}=[cellfiledata(ii).stimpath,cellfiledata(ii).stimfile];
   respfiles{ii}=[cellfiledata(ii).path,cellfiledata(ii).respfile];
end

fitfrac=0.15;
predfrac=0;

cellid=cellfiledata(1).cellid;
if cellfiledata(1).runclassid==5, % ie, WG
   respfmtcode=0;
   stimfmtcode=6;
   
   % resp load cmd. format: resploadcmd(respfile,p1,p2,...)
   resploadcmd='loadpypedata';
   resploadparms={'',1,0};
   
   % resp filter - eg, resample or smooth 
   respfiltercmd='';
   respfilterparms={};
   
   % stim load cmd. should be of format
   %   stimloadcmd(stimfile,startframe,endframe,p1,p2,...)
   stimloadcmd='loadwgstim';
   stimloadparms={4,4,5};
   
   % stim filter - eg, convert to phase-sep fourier domain.
   stimfiltercmd='';
   stimfilterparms={};
   
   kernfmt='wg';
   
   maxlag=[-8 15];
   resampfmt=1;     % ie bootstrapping or
   resampcount=20;   % number of resampled kernels
   
   dotSA=1;
   dosSA=1;  % 1 means first order decorrelation
   
   sffiltsigma=7;
   sfscount=30;
   neigs=0:4:(sfscount-1)*4;
   
   % need to decide:
   % which form of normalization? full, full regularized, ***1st order?***
   % what range of sigmas for threshold? * done? *
   % what rectification scheme in cellfit? *currently rect each chan*

else
   % set up for NR RC here....
   respfmtcode=0;
   
   % resp load cmd. format: resploadcmd(respfile,p1,p2,...)
   resploadcmd='loadpypedata';
   resploadparms={'',1,0};
   
   % resp filter - eg, resample or smooth 
   respfiltercmd='';
   respfilterparms={};
   
   % stim load cmd. should be of format
   %   stimloadcmd(stimfile,startframe,endframe,p1,p2,...)
   stimloadcmd='loadimfile';
   stimloadparms={};
   
   % stim filter - eg, convert to phase-sep fourier domain.
   stimfiltercmd='';
   stimfilterparms={};
   
   if strcmp(fmt,'sfft'),
      kernfmt='fft';
      stimfmtcode=7;
   else
      kernfmt='wav';
      stimfmtcode=5;
   end
   
   maxlag=[-8 15];
   resampfmt=1;     % ie bootstrapping or
   resampcount=20;   % number of resampled kernels
   
   dotSA=1;
   dosSA=3;  % 1 means first order decorrelation
             % 3 means preloaded (single sSA2?) 
   
   sffiltsigma=7;
   sfscount=30;
   neigs=0:2:(sfscount-1)*2;
   
end

% do the RC
cellxc;

% fit the kernels (best sfs and smoothing fudge factors to fit data)
cellfit;

% save important parms to movformat
movformat.iconside=iconside;
movformat.kernfmt=kernfmt;
movformat.daterun=date;
movformat.sigrange=sigrange;
movformat.sigfit=sigfit;
movformat.sfsidxtouse=sfsidxtouse;
movformat.sfsfit=sfsfit;
movformat.spacecount=spacecount;
movformat.maxlag=maxlag;
movformat.resampcount=resampcount;
movformat.stimfmtcode=stimfmtcode;
movformat.fitfrac=fitfrac;
movformat.predfrac=predfrac;
movformat.cellfiledata=cellfiledata;
movformat.boundary=boundary;
movformat.stimfiltercmd=stimfiltercmd;
movformat.stimfilterparms=stimfilterparms;
movformat.stimloadcmd=stimloadcmd;
movformat.stimloadparms=stimloadparms;
movformat.cellid=cellid;

% calculate step response
stephf=cumsum(hf,2);

steptime=sum(stephf,1);
maxidx=min(find(steptime==max(steptime)));
figure(2)
plot(steptime,'k-');
hold on
plot(maxidx,steptime(maxidx),'ro');
hold off
fprintf('peak of step response in bin %d. setting kern=stephf(:,%d)\n',...
        maxidx,maxidx);
kern=stephf(:,maxidx);

showkern(cat(3,hf(:,1:8),stephf(:,1:8)),kernfmt,iconside);
subplot(2,8,1);
title(sprintf('%s Impulse response',cellid));
subplot(2,8,9);
title(sprintf('%s Step response',cellid));


