% function [hf,stephf,movformat,kern]=rconline(rawids,fmt);
%
% rawids - vector of id tags from gDataRaw in cellDB. raw id is
%          listed next to the file name for each data set
% fmt - input transformation of stimulus to use for RC
%       for WG: use fmt='parm'
%       for NR: use fmt='sfft' for phase separated fourier domain
%               use fmt='wav' for wavelet domain
%
% returns:
% hf - STRF estimate in space X time dimensions
% stephf - step response of STRF (ie, cumsum(hf,2))
% movformat - parameters used for RC that can be used to guide
%             optinat stimulus generation
% kern - spatial only kernel. slice of stephf at maximum latency
%        (ie, greatest power before inhibition/adaptation kicks in)
%
% created SVD 5/13/02
%
function [hf,stephf,movformat,kern]=rconlineold(rawids,fmt);

% rawfileids gives entries into gDataRaw for data sets to analyze
if ~exist('rawids','var'),
   rawfileids=[596,597,598];
   %rawfileids=[592];
end
if ~exist('fmt','var'),
   fmt='pfft';
   fprintf('using fmt=''%s''\n',fmt);
end

% cellfileids give entries into sCellFile for data sets to analyze.
% if this is a new cell, dbraw2scellfile will load relevant info
% from gDataRaw and create appropriate entries in sCellFile in
% order to return valid indices
[cellfileids,cellfiledata]=dbraw2scellfile(rawids,fmt);

if length(cellfileids)==0,
   disp('no entries found in sCellFile!');
   return
end

dbopen;
global BATQUEUEID
BATQUEUEID=[];

%for ii=1:3,
%   figure(ii);
%end
drawnow;

if cellfiledata(1).runclassid==5, % ie, WG
   params.cellid=cellfiledata(1).cellid;
   params.stimfiles={};
   params.respfiles={};
   resplens=zeros(length(cellfiledata),1);
   for ii=1:length(cellfiledata),
      params.stimfiles{ii}=[cellfiledata(ii).stimpath,cellfiledata(ii).stimfile];
      params.respfiles{ii}=[cellfiledata(ii).path,cellfiledata(ii).respfile];
   end

   params.fitfrac=0.15;
   params.predfrac=0;
   
   params.respfmtcode=0;
   params.stimfmtcode=6;
   
   % resp load cmd. format: resploadcmd(respfile,p1,p2,...)
   params.resploadcmd='loadpypedata';
   params.resploadparms={'',1,0};
   
   % resp filter - eg, resample or smooth 
   params.respfiltercmd='';
   params.respfilterparms={};
   
   % stim load cmd. should be of format
   %   stimloadcmd(stimfile,startframe,endframe,p1,p2,...)
   params.stimloadcmd='loadwgstim';
   params.stimloadparms={4,4,5};
   params.stimiconside=[5 5 4 4];
   
   % stim filter - eg, convert to phase-sep fourier domain.
   params.stimfiltercmd='';
   params.stimfilterparms={};
   
   params.kernfmt='wg';
   
   params.minlag=-8;
   params.maxlag=15;
   params.resampfmt=1;     % ie bootstrapping or
   params.resampcount=20;   % number of resampled kernels
   
   params.dotSA=1;
   params.dosSA=1;  % 1 means first order decorrelation
   
   params.sffiltsigma=7;
   params.sfscount=30;
   params.sfsstep=4;
   %neigs=0:4:(sfscount-1)*4;
   
   % need to decide:
   % which form of normalization? full, full regularized, ***1st order?***
   % what range of sigmas for threshold? * done? *
   % what rectification scheme in cellfit? *currently rect each chan*
   
else
   batchid=35;
   
   fprintf('Loading parms for batchid=%d\n',batchid);
   
   sql=['SELECT * FROM sBatch WHERE id=',num2str(batchid)];
   batchdata=mysql(sql);
   params=batchdata;
   
   if strcmp(fmt,'sfft'),
      params.stimfmtcode=7;
      params.kernfmt='fft';
      stimiconside=[16 16];
      
   elseif strcmp(fmt,'downsamp'),
      params.stimfmtcode=2;
      params.kernfmt='space';
      stimiconside=[16 16];
      
   elseif strcmp(fmt,'wav'),
      params.stimfmtcode=5;
      params.kernfmt='wav';
      stimiconside=[5 5 4 4];
      
   elseif strcmp(fmt,'pfft'),
      params.stimfmtcode=8;
      params.kernfmt='pfft';
      stimiconside=[16 16];
      
   else
      disp('fmt must be: "space", "sfft", "pfft", or "wav" for NR data');
      return
   end
   
   % fill up params structure for passing to cellxcnodb
   params.cellid=cellfiledata(1).cellid;
   
   params.stimfiles={};
   params.respfiles={};
   resplens=zeros(length(cellfiledata),1);
   stimpix=zeros(length(cellfiledata),1);
   params.stimcrfs=zeros(length(cellfiledata),1);
   totlen=0;
   for ii=1:length(cellfiledata),
      params.stimfiles{ii}=[cellfiledata(ii).stimpath,...
                    cellfiledata(ii).stimfile];
      params.respfiles{ii}=[cellfiledata(ii).path,...
                    cellfiledata(ii).respfile];
      resplens(ii)=cellfiledata(ii).resplen;
      tstimpix=strsep(cellfiledata(ii).stimiconside,',');
      if length(tstimpix)>0,
         stimpix(ii)=tstimpix{1};
      end
      if cellfiledata(ii).stimfilecrf>0,
         params.stimcrfs(ii)=cellfiledata(ii).stimfilecrf;
      else
         sql=['SELECT * FROM gCellMaster WHERE cellid="',params.cellid,'"'];
         celldata=mysql(sql);
         params.stimcrfs(ii)=stimpix(ii)./celldata.rfsize;
      end
   end
   
   params.batch=batchdata.id;
   params.stimwindowsize=cellfiledata(1).stimwindowsize;
   
   % these entries in batchdata need to be parsed
   params.resploadparms=strsep(batchdata.resploadparms,',');
   params.respfilterparms=strsep(batchdata.respfilterparms,',');   
   params.stimloadparms=strsep(batchdata.stimloadparms,',');
   params.stimfilterparms=strsep(batchdata.stimfilterparms,',');
   
   params.docellfit2=0;
   params.shrinkage=1;   % rather than just thresholding
   params.repexclude=0;
   
   params.outfile='/tmp/rconline.mat';
   params.zipoutfile=0;
   
end

cellxcnodb(params);

load(params.outfile);
movformat=params;
movformat.cellfiledata=cellfiledata;

nlidx=1;
hf=strf(nlidx).h;
xc=xc(:,:,nlidx,1);
predxc=predxc(:,nlidx,1);
sigfit=strf(nlidx).parms.sigfit;
sfsfit=strf(nlidx).parms.sfsfit;

stephf=cumsum(hf,2);
steptime=sum(stephf,1);
maxidx=min(find(steptime==max(steptime)));
kern=stephf(:,maxidx);

fprintf('peak of step response in bin %d. (nlidx=%d) kern=stephf(:,%d)\n',...
        maxidx,nlidx,maxidx);

save /tmp/rconline.mat hf stephf movformat kern
fprintf('hf stephf movformat kern saved to /tmp/rconline.mat\n');


