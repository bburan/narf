function viewpfth(cellid,batch)

[cellfiledata,times,batchdata]=cellfiletimes(cellid,batch);
starttimes=times(1).start;
stoptimes=times(1).stop;
fitfile=times(2).fileidx;
fitstartframe=times(2).start;
fitstopframe=times(2).stop;
predfile=times(3).fileidx;
predstartframe=times(3).start;
predstopframe=times(3).stop;

stimfiles={};
respfiles={};
resplens=zeros(length(cellfiledata),1);
stimpix=zeros(length(cellfiledata),1);
stimcrfs=zeros(length(cellfiledata),1);
totlen=0;
for ii=1:length(cellfiledata),
   stimfiles{ii}=[cellfiledata(ii).stimpath,cellfiledata(ii).stimfile];
   respfiles{ii}=[cellfiledata(ii).path,cellfiledata(ii).respfile];
   resplens(ii)=cellfiledata(ii).resplen;
   tstimpix=strsep(cellfiledata(ii).stimiconside,',');
   stimpix(ii)=tstimpix{1};
   stimcrfs(ii)=cellfiledata(ii).stimfilecrf;
end

stimfmtcode=batchdata.stimfmtcode;
respfmtcode=batchdata.respfmtcode;

% resp load cmd. format: resploadcmd(respfile,p1,p2,...)
resploadcmd=batchdata.resploadcmd;
resploadparms=strsep(batchdata.resploadparms,',');

% resp filter - eg, resample or smooth 
respfiltercmd=batchdata.respfiltercmd;
respfilterparms=strsep(batchdata.respfilterparms,',');

% stim load cmd. should be of format
%   stimloadcmd(stimfile,startframe,endframe,p1,p2,...)
stimloadcmd=batchdata.stimloadcmd;
stimloadparms=strsep(batchdata.stimloadparms,',');

% stim filter - eg, convert to phase-sep fourier domain.
stimfiltercmd=batchdata.stimfiltercmd;
stimfilterparms=strsep(batchdata.stimfilterparms,',');
stimwindowcrf=batchdata.stimwindowcrf;

kernfmt=batchdata.kernfmt;

maxlag=[batchdata.minlag batchdata.maxlag];
resampfmt=batchdata.resampfmt;     % ie bootstrapping or
resampcount=batchdata.resampcount; % number of resampled kernels
                                   % something else
dotSA=batchdata.decorrtime;
dosSA=batchdata.decorrspace;

neigs=0:batchdata.sfsstep:(batchdata.sfsstep*batchdata.sfscount-1);
sfscount=batchdata.sfscount;
sffiltsigma=batchdata.sffiltsigma;
nloutparm=batchdata.nloutparm;

fidx=1;

% check to see if stim and response exist
if ~checkbic & ~exist(respfiles{fidx},'file'),
   disp([respfiles{fidx},' missing!!']);
   return
end
if ~checkbic & ~exist(stimfiles{fidx},'file'),
   disp([stimfiles{fidx},' missing!!']);
   return
end

%
% load the response data
%
disp(sprintf('Loading response: %s...',respfiles{fidx}));
% resp=resploadatt(respfile,respvarname);
% resp is time X space(phase) X attcode
%bresp=feval(resploadcmd,respfiles{fidx},resploadparms{:});   
bresp=resploadatt(respfiles{fidx},'pfth',1,15,1);   
cresp=resploadatt(respfiles{fidx},'pfth',1,6,1);   

% filter response (eg resample, pick attentional state, etc)
if ~isempty(respfiltercmd),
   %bresp=feval(respfiltercmd,bresp,respfilterparms{:});
end
rsize=size(bresp);

bstim=feval(stimloadcmd,stimfiles{fidx},1,0,stimloadparms{:});
   
% filter stimulus segment if selected
%if ~isempty(stimfiltercmd),
%   bstim=feval(stimfiltercmd,bstim,stimfilterparms{:});
%end
% reshape to space X time if necessary
iconside=size(bstim);
iconside=iconside(1:(end-1));
%f length(iconside)>=2,
%  % reshape all spatial dims into one
%  bstim=reshape(bstim,prod(iconside),size(bstim,length(iconside)+1));
%nd
%stim=bstim'; % take transpose to put time in rows
spacecount=prod(iconside);

keyboard

if 1| ~ismember(stimfmtcode,[4]), % ie, not mfilt
   [celldata,rawdata]=getfreefiles(cellid);
   rawdata=rawdata(1);
   
   % figure out window parameters used for stimulus generation
   movformat.stimwindowsize=celldata.rfsize;
   if movformat.stimwindowsize>200,
      movformat.stimwindowsize=200;
   end
   movformat.offx=celldata.xoffset;
   movformat.offy=celldata.yoffset;
   movformat.rfsize=celldata.rfsize;
   
   freefile=[rawdata.resppath,rawdata.matlabfile];
   movformat.rawid=rawdata.id;
   
   z=zload([freefile,'.gz'],checkbic);
   if isempty(z),
      return
   end
   
   % maintain backward compatibility
   if isfield(z,'data'),
      data=z.data;
   elseif isfield(z,'z'),
      data=z.z;
   end
   clear z
   
   % figure out target ids. need to be certain that this matches respfile!
   if isfield(data,'saclist'),
      targlist=unique(data.saclist.sample)+1;
   else
      tlist=[];
      for ii=1:length(data.trials),
         if data.trials(ii).valid,
            tlist=[tlist data.trials(ii).f(1,3)];
         end
      end
      targlist=unique(tlist);
      
      saclist=[];
   end
   
   fprintf('Targets id''d:');
   for ii=1:length(targlist),
      fprintf(' %d',targlist(ii));
   end
   fprintf('\n');
   
   disp('generating patches for matched filter analysis...');
   
   if stimfmtcode==5 | stimfmtcode==4,
      iconside=[movformat.stimwindowsize movformat.stimwindowsize];
   else
      iconside=strsep(cellfiledata(1).stimiconside,',');
      iconside=cat(2,iconside{:});
   end
   
   imagepix=size(data.images{1},1);
   patchcount=length(data.images);
   patchlist=1:patchcount;
   if isfield(data.trials(1),'params') & ...
         isfield(data.trials(1).params,'smooth1'),
      alphamask=makealphamask(imagepix,data.trials(1).params.smooth1,...
                              data.trials(1).params.smooth2);
      patchpix=data.trials(1).params.smooth2*2;
      targpatches=ones([iconside(:)',patchcount]) .* alphamask(1,1);
   else
      alphamask=ones(imagepix);
      targpatches=ones([iconside(:)',patchcount]).*data.images{1}(1,1);
   end
   patchpix=movformat.stimwindowsize;
   bigpatches=zeros(patchpix,patchpix,patchcount);
   
   for ii=1:patchcount,
      
      % insert patch in window sized to rf.
      if imagepix>=patchpix,
         m=floor((imagepix+1)/2);
         timage=alphamask .* (data.images{patchlist(ii)}/128-1);
         timage=timage(m-patchpix/2+1:m+patchpix/2,...
                       m-patchpix/2+1:m+patchpix/2);
      else
         m=floor((patchpix+1)/2);
         timage=zeros(patchpix);
         timage(m-imagepix/2+1:m+imagepix/2,...
                m-imagepix/2+1:m+imagepix/2)=...
            alphamask .* data.images{patchlist(ii)}/128-1;
      end
      
      fmask=ones(patchpix);
      croprange=[1 1 patchpix patchpix];
      targpatches(:,:,ii)=movresize(timage,iconside(1),fmask,croprange,0,1);
      bigpatches(:,:,ii)=movresize(timage,patchpix,fmask,croprange,0,1);
   end
end

brmean=nanmean(bresp(:,:,1));
crmean=nanmean(cresp(:,:,1));
bsmean=mean(bstim,1);

brmax=min(find(brmean==max(brmean)));
fprintf('brmax=%d\n',brmax);

attcount=size(bresp,3);
pcount=10;
rlen=rsize(1);

figure
clf
for attidx=1:attcount,
   [bs,bsidx]=sort(bresp(:,brmax,attidx));
   rlen=max(find(~isnan(bs)));
   for pidx=1:pcount,
      if pidx>1,
         imidx=rlen-pcount+pidx;
         im=bstim(:,:,bsidx(imidx));
      elseif attidx>1,
         %im=bigpatches(:,:,targlist(attidx-1));
         im=targpatches(:,:,targlist(attidx-1));
      else
         im=zeros(iconside);
      end
      subplot(attcount,pcount,(attidx-1)*pcount+pidx);
      imagesc(im,[-1 1]);
      axis off
      axis image
      if (pidx==1),
         title(sprintf('%s attidx=%d',cellid,attidx));
      else
         title(sprintf('r=%.3f',bs(imidx)));
      end
   end
end
colormap(gray);
