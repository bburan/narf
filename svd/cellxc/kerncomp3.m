% kerncomp.m : compare kernels from different attentional states to
% see how similar they are

disp('kerncomp3.m');

% parms for XC / decorr
if ~exist('maxlag'),
   %maxlag=[-10 15];
   maxlag=[0 0];
end
if ~exist('meansub'),
   meansub=1;
end
if ~exist('boundary'),
   boundary='zero';
end
if ~exist('mineigs'),
   mineigs=[];
end
% movstep: size of movie segments to send to movxc -- want this to
% be as large as possible without straining computer's memory
if ~exist('movstep'),
   movstep=50000;
end

filecount=length(respfiles);

% figure out separate start and stop times for each attentional state!
rvalid={};

if times(1).stop(1)-times(1).start(1) > 40000,
   PRELOADMOV=2;
else
   PRELOADMOV=1;
end

DOREGNORM=0;

spacelims=20;
spacelimp=20;
PFILT=0.010;  % p value (any attention) to be considered stimulus-modulated
PATT=0.008;  % p<0.008 any pair has p<0.05 attention modulation

% noisecount already defined... number of fake "attentional" states
% attcount: hacky parameter telling how many attentional conditions
% there are (4 targets plus "all/-1"). adjusted 
attcount=5;

%
% FIGURE OUT WHAT THIS HAS TO TO WITH THE TARGETS!
%

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

fprintf('loading patch data from %s...\n',freefile);
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

%keyboard

for ii=1:patchcount,
   
   % insert patch in window sized to rf.
   if imagepix>=patchpix,
      m=floor((imagepix+1)/2);
      timage=alphamask .* (data.images{patchlist(ii)}) + ...
             (1-alphamask) .* data.bgpix;
      timage=timage(m-patchpix/2+1:m+patchpix/2,...
                    m-patchpix/2+1:m+patchpix/2);
   else
      m=floor((patchpix+1)/2);
      timage=ones(patchpix) .* data.images{patchlist(ii)}(1);
      timage(m-imagepix/2+1:m+imagepix/2,...
             m-imagepix/2+1:m+imagepix/2)=...
         alphamask .* data.images{patchlist(ii)} + ...
         (1-alphamask) .* data.bgpix;
   end
   
   fmask=ones(patchpix);
   croprange=[1 1 patchpix patchpix];
   targpatches(:,:,ii)=movresize(timage,iconside(1),fmask,croprange,0,1);
   bigpatches(:,:,ii)=movresize(timage,patchpix,fmask,croprange,0,1);
end

% simulate what happened to stimulus during RC:

% does it need to be resized again???
if stimfmtcode==5,
   wavside=strsep(cellfiledata(1).stimiconside,',');
   wavside=cat(2,wavside{:});
   spacebincount=wavside(1);
   obincount=wavside(3);
   sbincount=wavside(4);
   targpatches=mov2wav(targpatches,spacebincount,obincount,sbincount);
   targpatches=reshape(targpatches,prod(wavside),patchcount);
elseif stimfmtcode==4,
   targpatches=eye(patchcount);
elseif length(stimloadparms)>0,
   % hacked from loadimfile...
   [fmask,crop]=movfmask(iconside(1),stimloadparms{2}./stimloadparms{1},...
                         stimloadparms{3}*iconside(1)./stimloadparms{1});
   
   targpatches=movresize(targpatches,stimloadparms{3},fmask,crop,0,1);
end

% and do the appropriate stim filter:
if ~isempty(stimfiltercmd),
   fpatches=feval(stimfiltercmd,targpatches,stimfilterparms{:});
else
   fpatches=targpatches;
end
iconside=size(fpatches);
iconside=iconside(1:(end-1));
if length(iconside)>=2,
   % reshape all spatial dims into one
   fpatches=reshape(fpatches,prod(iconside),...
                    size(fpatches,length(iconside)+1));
end
mpatches=[mean(fpatches,2),fpatches(:,targlist)];

% clear big variables to save memory
clear data z

% for each file:
for fidx=1:filecount,
   
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
   bresp=feval(resploadcmd,respfiles{fidx},resploadparms{:});   
   
   % filter response (eg resample, pick attentional state, etc)
   if ~isempty(respfiltercmd),
      bresp=feval(respfiltercmd,bresp,respfilterparms{:});
   end
   
   bs=size(bresp);
   attcount=bs(3);
   
   % figure out size and mean of response
   rsize=size(bresp);
   blen=rsize(1);
   respcount=rsize(2);
   brmean=zeros(1,respcount,attcount);
   for respidx=1:respcount,
      brmean(:,respidx,:)=nanmean(bresp(:,respidx,:));
   end
   
   % adjust attcount to actual number of attention states
   if fidx==1,
      attcount=size(bresp,3);  % dim 3 is number attentional states.
      resplens=zeros(filecount,attcount);
      starttimes=zeros(filecount,attcount);
      stoptimes=zeros(filecount,attcount);
   end
   for attidx=1:attcount,
      rvalid{fidx,attidx}=find(~isnan(squeeze(bresp(:,1,attidx))));
      resplens(fidx,attidx)=length(rvalid{fidx,attidx});
      if resplens(fidx,attidx)>0, 
         starttimes(fidx,attidx)=rvalid{fidx,attidx}(1);
         stoptimes(fidx,attidx)=rvalid{fidx,attidx}(resplens(fidx,attidx));
      end
   end
   
   if PRELOADMOV==1,
      bstim=feval(stimloadcmd,stimfiles{fidx},1,0,stimloadparms{:});
      
      % filter stimulus segment if selected
      if ~isempty(stimfiltercmd),
         bstim=feval(stimfiltercmd,bstim,stimfilterparms{:});
      end
      % reshape to space X time if necessary
      iconside=size(bstim);
      iconside=iconside(1:(end-1));
      if length(iconside)>=2,
         % reshape all spatial dims into one
         bstim=reshape(bstim,prod(iconside),size(bstim,length(iconside)+1));
      end
      bstim=bstim'; % take transpose to put time in rows
      spacecount=size(bstim,2);
      bsmean=mean(bstim,1);
      
      tSA1=zeros(diff(maxlag)*2+1,1);
      for xx=1:spacecount,
         if meansub,
            tstim=bstim(:,xx)-mean(bstim(:,xx));
         else
            tstim=bstim(:,xx);
         end
         tSA1=tSA1+xcorr(tstim,diff(maxlag),'biased')./spacecount;
      end
      U=diff(maxlag)+1;
      for resampidx=1:resampcount,
         ttSA=zeros(U);
         for uu=1:U,
            ttSA(:,uu)=tSA1((U+1-uu):(U*2-uu));
         end
         ttSAinv=svdinv(ttSA,0.0001);
      end
   elseif PRELOADMOV==2 & ~checkbic,
      disp('Copying stimulus file to local drive');
      tstimfile=[getenv('TMP'),'/',basename(stimfiles{fidx})];
      unix(['cp ',stimfiles{fidx},' ',tstimfile]);
   else
      tstimfile=stimfiles{fidx};
   end
   
   bsmean=mean(bstim,1);
   bstim=bstim-repmat(bsmean,[blen,1]);
   
   %
   % BEGIN STRF SECTION
   %
   % 0. compute STRF in the PC domain
   disp('0. Testing for linear modulation...');
   
   %keyboard
   
   segcount=20;
   SR=zeros(spacecount,respcount,attcount,noisecount+1);
   Hpc=zeros(spacelims,respcount,attcount,noisecount+1);
   sSA=zeros(spacecount,spacecount,attcount);
   
   rmean=reshape(nanmean(bresp),respcount,attcount);
   rstd=zeros(respcount,attcount);
   for attidx=1:attcount,
      rstd(:,attidx)=nanstd(bresp(:,:,attidx))';
   end
   pstrf=ones(spacelims,respcount,attcount);
   for attidx=1:attcount,
      
      % compute spatial correlation for stim in attidx
      if attidx==1,
         aokidx=find(sum(~isnan(bresp(:,end,2:end)),3));
      else
         aokidx=find(~isnan(bresp(:,end,attidx)));
      end
      stim=bstim(aokidx,:);
      smean=mean(stim,1);
      stim=stim-repmat(smean,[size(stim,1),1]);
      
      sSA(:,:,attidx)=stim'*stim ./ length(aokidx);
      
      % compute matrix to convert everyone into same, decorrelated space
      if attidx==1,
         [ua,sa,va]=svd(sSA(:,:,attidx));
         B=diag(1./diag(sa))*ua';
         B=B(1:spacelims,:);
      else
         % transform strf from attention state idx into PC domain
         % of all att, decorrelating as needed.
         [ut,st,vt]=svd(sSA(:,:,attidx));
         
         % ssa=u*s*v'
         % ssa_pc=ua'* (u * s * v') * ua
         % ssa^-1=v*s^-1*u'
         % ssa_pc^-1 = ua'* v* s^-1 * u'* ua
         
         sst=1./diag(st);
         %sst(end)=0;
         B=ua' * vt * diag(sst) * ut';
         B=B(1:spacelims,:);
      end
      
      for respidx=1:respcount,
         fprintf('attidx=%d respidx=%d...\n',attidx,respidx);
         okidx=find(~isnan(bresp(:,respidx,attidx)));
         if ~isempty(okidx),
            r=bresp(okidx,respidx,attidx)-rmean(respidx,attidx);
            stim=bstim(okidx,:);
            slen=size(stim,1);
            smean=mean(stim,1);
            stim=stim-repmat(smean,[slen,1]);
            
            SR(:,respidx,attidx,1)=stim'*r./slen;
            Hpc(:,respidx,attidx,1)=B*SR(:,respidx,attidx,1);
            
            rnoise=zeros(slen,noisecount);
            for noiseidx=1:noisecount,
               [x,setidx]=sort(rand(size(r)));
               rnoise(:,noiseidx)=r(setidx);
            end
            tSR=stim'*rnoise./slen;
            tHpc=B*tSR;
            
            SR(:,respidx,attidx,2:end)=reshape(tSR,spacecount,1,1,noisecount);
            Hpc(:,respidx,attidx,2:end)=reshape(tHpc,spacelims,1,1,noisecount);
            
            for eigidx=1:spacelims,
               nn=length(find(abs(Hpc(eigidx,respidx,attidx,2:end))>=...
                              abs(Hpc(eigidx,respidx,attidx,1))));
               
               pstrf(eigidx,respidx,attidx)=(nn+1)./(noisecount+1);
            end
         end
         if exist('BATQUEUEID') & BATQUEUEID>0,
            dbsetqueue(BATQUEUEID,respidx+attidx*10);
         end
      end
   end
   
   % record original number of spatial dims tested
   startspacelims=spacelims;
   
   for eigidx=1:spacelims,  % spacecount
      if Hpc(eigidx,end,1,1)<0,  % sum(ua(:,eigidx))<0,
         ua(:,eigidx)=-ua(:,eigidx);
         va(:,eigidx)=-va(:,eigidx);
         Hpc(eigidx,:,:,:)=-Hpc(eigidx,:,:,:);
      end
   end
   
   % save original STRF
   Hpcfull=Hpc;
   %Hpc=Hpcfull;
   
   if 1, 
      
      % rank strf coefficients by snr.
      Hm=Hpc(:,:,:,1);
      Hs=mean(std(Hpc(:,:,2:end,2:end),1,4),3);
      
      resprange=1:respcount; 
      attrange=1;
      snr=max(max(abs(Hm(:,resprange,attrange))./Hs(:,resprange,attrange),...
                  [],2),[],3);
      [ssnr,spaceuses]=sort(-snr);
      ssnr=-ssnr;
      
      % reorder by snr
      Hpc=Hpc(spaceuses,:,:,:);
      Hs=Hs(spaceuses,:,:);
      
      % factor it apart.  scale each attentional kernel by snr of
      % all-att state.
      respidx=1;
      H=squeeze(Hpc(:,respidx,2:end,1))./repmat(Hs(:,1),[1 attcount-1]);
      
      [h0,s0,att0]=svd(H);
      ss=diag(s0);
      for ii=1:size(att0,2),
         if sum(att0(:,ii))<0,
            att0(:,ii)=-att0(:,ii);
            h0(:,ii)=-h0(:,ii);
         end
      end
      
      %keyboard
      
   else
      % find valid spatial dims in STRF
      ptest=permute(pstrf,[3 1 2]);
      if min(ptest(:))>PFILT,
         PFILT=min(ptest(:))+PFILT;
      end
      spaceuses=find(sum(sum(ptest<PFILT,1),3)>0);
      spacelims=length(spaceuses);
      
      fprintf('%d significant dims: ',spacelims);
      for eigidx=spaceuses(:)',
         fprintf('%d ',eigidx);
      end
      fprintf(' (p<%.3f)\n',PFILT);
                   
      % remove unmodulated dimensions from the mix
      Hpc=Hpc(spaceuses,:,:,:);
   end
   
   
   % 1. determine significantly attention-modulated spatial channels
   disp('1. Calculating STRF attention-modulation T2 and significance...');
   
   % is attention comparison in/out or pairwise in/in? currently
   % pairwise seems best
   DOPAIRS=1;
   if DOPAIRS,
      attpairs=[1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
      dopairs=find(attpairs(:,1)+1<=attcount & attpairs(:,2)+1<=attcount);
      attpairs=attpairs(dopairs,:);
      attdocount=size(attpairs,1);
   else
      attdocount=4;
   end
   
   T2strf=zeros(attdocount,spacelims,respcount,noisecount+1,2);
   pstrfatt=ones(attdocount,spacelims,respcount);
   ppcatt=ones(attdocount,spacelims,respcount);
   for pidx=1:attdocount,
      a1=attpairs(pidx,1)+1;
      a2=attpairs(pidx,2)+1;
      fprintf('paircomp %d v %d',a1,a2);
      for respidx=1:respcount,
         fprintf(' r=%d',respidx);
         
         % then estimate noise kernels and calculate T2 for each
         % of them..??..?? ugly, slow, but reasonable?
         a1idx=find(~isnan(bresp(:,respidx,a1)));
         a2idx=find(~isnan(bresp(:,respidx,a2)));
         aidx=[a1idx;a2idx];
         r=[bresp(a1idx,respidx,a1); bresp(a2idx,respidx,a2)];
         
         nSR=zeros(spacecount,noisecount+1,2);
         nHpc=zeros(spacelims,noisecount+1,2);
         nSR(:,1,:)=SR(:,respidx,[a1 a2],1);
         nHpc(:,1,:)=Hpc(:,respidx,[a1 a2],1);
         
         noisecol=zeros(length(aidx),noisecount);
         rcol1=zeros(length(a1idx),noisecount);
         rcol2=zeros(length(a2idx),noisecount);
         for noiseidx=1:noisecount,
            [x,setidx]=sort(rand(size(aidx)));
            n1idx=setidx(1:length(a1idx));
            n2idx=setidx(length(a1idx)+1:end);
            
            noisecol(n1idx,noiseidx)=1;
            rcol1(:,noiseidx)=r(n1idx);
            rcol2(:,noiseidx)=r(n2idx);
            
            rnoise=r(n1idx)-mean(r(n1idx));
            stim=bstim(aidx(n1idx),:);
            slen=size(stim,1);
            smean=mean(stim,1);
            stim=stim-repmat(smean,[slen,1]);
            
            nSR(:,noiseidx+1,1)=stim'*rnoise./slen;
            nSA=stim'*stim ./ slen;
            [ut,st,vt]=svd(nSA);
            sst=1./diag(st);
            %sst(end)=0;
            B=ua(:,spaceuses)' * vt * diag(sst) * ut';
            nHpc(:,noiseidx+1,1)=B*nSR(:,noiseidx+1,1);
            
            rnoise=r(n2idx)-mean(r(n2idx));
            stim=bstim(aidx(n2idx),:);
            slen=size(stim,1);
            smean=mean(stim,1);
            stim=stim-repmat(smean,[slen,1]);
            
            nSR(:,noiseidx+1,2)=stim'*rnoise./slen;
            nSA=stim'*stim ./ slen;
            [ut,st,vt]=svd(nSA);
            sst=1./diag(st);
            %sst(end)=0;
            B=ua(:,spaceuses)' * vt * diag(sst) * ut';
            nHpc(:,noiseidx+1,2)=B*nSR(:,noiseidx+1,2);
            
            if exist('BATQUEUEID') & BATQUEUEID>0,
               dbsetqueue(BATQUEUEID,1000*pidx+500*respidx+noiseidx);
            end
         end
         
         if 1
         % calc T2 for actual from Hpc
         % with error bars 2:end actual attentional conditions
         for eigidx=1:spacelims,
            
            % hackish? borrow error bars (xstd) from actual att
            % kernels for all the noise-attention kernels
            
            x1std=std(Hpc(1:eigidx,respidx,a1,2:end),0,4);
            x2std=std(Hpc(1:eigidx,respidx,a2,2:end),0,4);
            xstd=(x1std+x2std)./2;
            xstd(find(xstd==0))=1;
            
            for noiseidx=1:noisecount+1,
               x1=nHpc(1:eigidx,noiseidx,1);
               x2=nHpc(1:eigidx,noiseidx,2);
               
               T2strf(pidx,eigidx,respidx,noiseidx,1)=...
                   sum(abs(x1(:)-x2(:))./xstd(:),1);
               
               T2strf(pidx,eigidx,respidx,noiseidx,2)=...
                   abs(x1(eigidx)-x2(eigidx))./xstd(eigidx);
            end
         end
         end
         
         if 0,
         for eigidx=1:spacelims,
            for noiseidx=1:noisecount+1,
               % for each pair, exclude 1, compute T2 for reama
               x1=nHpc(1:eigidx,noiseidx,1);
               x2=nHpc(1:eigidx,noiseidx,2);
               
               E=x1-x2;
               
               X=nHpc(1:eigidx,[1:noiseidx-1 noiseidx+1:end],:);
               mX=mean(X,2);
               X=X-repmat(mX,[1 noisecount 1]);
               V=(X(:,:,1)*X(:,:,1)' + X(:,:,2)*X(:,:,2)') ./ ...
                 (noisecount-2);
               
               % if V is nice and nonsingular...
               Vinv=V^-1;
               
               % pseudo inverse... seems to flail in extremes
               %Vinv=svdinv(V,0.000001);
               
               % ignore correlations between channels?
               % amplifies noise but may be the only stable option
               %sV=1./diag(V);
               %sV(find(sV>median(sV)*10^6))=0;
               %Vinv=diag(sV);
               
               T2strf(pidx,eigidx,respidx,noiseidx)=...
                   E' * Vinv * E ./ size(X,1);
            end
         end
         end
         
         for eigidx=1:spacelims
            % compute pstrfatt value for current pair,eigextent,resp
            nn=length(find(abs(T2strf(pidx,eigidx,respidx,2:end))>=...
                           abs(T2strf(pidx,eigidx,respidx,1))));
            
            pstrfatt(pidx,eigidx,respidx)=(nn+1)./(noisecount+1);
         end
         
         if exist('BATQUEUEID') & BATQUEUEID>0,
            dbsetqueue(BATQUEUEID,200+10*pidx+respidx);
         end
      end
      fprintf('\n');
   end
   
   respidx=min([respcount,7]);
   
   figure(1);
   clf
   ps={'k-','r--','b--','g--','c--'};
   ps2={'kx','ro','bo','go','co'};
   labels={'A','1','2','3','4'};
   ll=zeros(attcount,1);
   
   subplot(2,1,1);
   hold on
   attshowset=1:attcount;
   for attidx=attshowset,
      a=errorbar(Hpcfull(:,respidx,attidx,1),...
               std(Hpcfull(:,respidx,attidx,2:end),0,4), ...
               ps{attidx});
      ll(attidx)=a(1);
      plot(spaceuses,Hpc(:,respidx,attidx,1),ps2{attidx});
   end
   aa=axis;
   axis([0 size(Hpcfull,1)+1 aa(3) aa(4)]);
   hold off
   title(sprintf('%s PC domain STRF',cellid));
   xlabel('pc idx');
   legend(ll(find(ll)),labels{find(ll)});
   
   hs=subplot(2,1,2);
   %imagesc(spaceuses,1:attdocount,pstrfatt(:,:,respidx),[0 1])
   imagesc(pstrfatt(:,:,respidx),[0 1]);
   axis xy
   hold on
   for pidx=1:attdocount,
      lidx=find(pstrfatt(pidx,:,respidx)<PATT);
      if length(lidx)>0,
         plot(lidx,pidx,'kx');
      end
   end
   hold off
   ylabel('1-2  1-3  1-4  2-3  2-4  3-4');
   xlabel(sprintf('signif pc idx (p<%.3f)',PATT));
   set(hs,'YTickLabel',[]);
   ttick=floor(str2num(get(hs,'XTickLabel')));
   ttick(find(ttick==0))=1;
   set(hs,'XTickLabel',spaceuses(ttick));
   colormap(flipud(hot));
   colorbar
   drawnow
      
   %
   % BEGIN PRF SECTION
   %
   
   [ua,sa,va]=svd(sSA(:,:,1));
   for eigidx=1:spacelimp,
      if sum(ua(:,eigidx))<0,
         ua(:,eigidx)=-ua(:,eigidx);
         va(:,eigidx)=-va(:,eigidx);
      end
   end
   
   % convert stim to eig domain
   if 1,
      % COMPUTE ERF IN REDUCED PC-DOMAIN SET
      eigstim=bstim * ua(:,1:spacelimp);
      eigmpatches=(mpatches'-repmat(bsmean,[attcount,1]))*ua(:,1:spacelimp);
      eigfpatches=(fpatches'-repmat(bsmean,[patchcount,1]))*ua(:,1:spacelimp);
   else
      % COMPUTE ERF IN "RAW" spatial dimensions
      spacelimp=128;
      eigstim=bstim;
      
      eigmpatches=(mpatches'-repmat(bsmean,[attcount,1]));
      eigfpatches=(fpatches'-repmat(bsmean,[patchcount,1]));
      disp('BIG PRF');
   end
   
   bincount=7;
   rbincount=zeros(bincount,respcount,attcount,spacelimp);
   mpbinned=zeros(attcount,spacelimp);
   sbin=zeros(bincount,spacelimp);
   binstim=zeros(size(eigstim));
   
   % cycle through each stimulus eigenvector
   %rbinset=cell(bincount,respcount,attcount,spacelimp);
   for eigidx=1:spacelimp ,
      fprintf('binning r: eigidx=%d ',eigidx);
      
      % coefficients of the stimulus for eigenvector eigidx
      sbinned=eigstim(:,eigidx);
      [ss,si]=sort(sbinned);
      
      tr=bresp .* repmat(eigstim(:,eigidx),[1 respcount,attcount]);
      for bb=1:bincount,
         
         % determine mean value of the bb-th fraction of the
         % stimulus coeffiecients
         sbin(bb,eigidx)=mean(sbinned(si(round((bb-1)*blen/bincount+1):...
                                           round(bb*blen/bincount))));
         
         % figure out which fixations are associated with that
         % stimulus bin
         sbinned(si(round((bb-1)*blen/bincount+1):round(bb*blen/bincount)))=bb;
         sbbidx=find(sbinned==bb);
         
         % determine number of responses falling in the stimulus
         % bin for each latency and attentional state
         rbincount(bb,:,:,eigidx)=sum(~isnan(tr(sbbidx,:,:)),1);
         %for respidx=1:respcount,
         %   for attidx=1:attcount,
         %      sbbokidx=sbbidx(find(~isnan(tr(sbbidx,respidx,attidx))));
         %      rbinset{bb,respidx,attidx,eigidx}=...
         %          tr(sbbokidx,respidx,attidx)./sbin(bb,eigidx);
         %   end
         %end
         
         fprintf('.');
         if exist('BATQUEUEID') & BATQUEUEID>0,
            dbsetqueue(BATQUEUEID,eigidx*100+bb);
         end
      end
      
      binstim(:,eigidx)=sbinned;
      
      % assign each eigenvector coefficient of the search targets
      % to a spatial bin
      for attidx=1:attcount,
         tt=abs(eigmpatches(attidx,eigidx)-sbin(:,eigidx));
         mpbinned(attidx,eigidx)=find(tt==min(tt));
      end
      
      fprintf('\n');
   end
   
   rmean=reshape(nanmean(bresp),respcount,attcount);
   rstd=zeros(respcount,attcount);
   for attidx=1:attcount,
      rstd(:,attidx)=nanstd(bresp(:,:,attidx))';
   end
   
   % 2. determine significantly modulated spatial channels
   % (independent of attention)
   %
   % compute chi^2 for actual erfs and noisecount sets randomized
   % over stimulus bin
   
   disp('2. Testing for attention-independent modulation in PRF...');
   % record mean and std of each bin for actual kernels
   rbinned=zeros(bincount,respcount,attcount,spacelimp);
   rbinnedstd=zeros(bincount,respcount,attcount,spacelimp);
   rm=zeros(bincount,1);
   rs=zeros(bincount,1);
   
   % chi^2 records for each actual and noise kernel, perf measures
   % p value for each (latidx, attidx, eigidx)
   fchi2=zeros(attcount,spacelimp,noisecount+1,respcount);
   perf=ones(attcount,spacelimp,respcount);
   
   for respidx=1:respcount,
      fprintf('respidx=%d\n',respidx);
      for attidx=1:attcount,
if 0
         % THIS IS LOESS EXPERIEMENTAL STUFF
         rvalididx=find(~isnan(bresp(:,respidx,attidx)));
         nerf=zeros(bincount,spacelimp,noisecount+1);
         nerfstd=zeros(bincount,spacelimp,noisecount+1);
         tbinstim=binstim(rvalididx,:);
         teigstim=eigstim(rvalididx,:);
         tr=bresp(rvalididx,respidx,attidx);
         
         %g = loess(x,y,newx,alpha,lambda,robustFlag)
         loalpha=0.6;
         lolambda=1;
         resampcount=20;
         slen=length(rvalididx);
         rstart=floor(linspace(0,slen,resampcount+1));
         
         for eigidx=1:spacelimp,
            terf=zeros(bincount,resampcount);
            for resampidx=1:resampcount,
               %trange=[1:rstart(resampidx) (rstart(resampidx+1)+1):slen];
               trange=[rstart(resampidx)+1:rstart(resampidx+1)];
               terf(:,resampidx)=loess(teigstim(trange,eigidx),tr(trange),...
                                       sbin(:,eigidx),loalpha,lolambda)';
            end
            rbinned(:,respidx,attidx,eigidx)=mean(terf,2);
            %rbinnedstd(:,respidx,attidx,eigidx)=std(terf,0,2);
            rbinnedstd(:,respidx,attidx,eigidx)= ...
                std(terf,0,2)./sqrt(resampcount);
            tt{attidx}=terf;
         end
         
         
         
         nerf(:,:,1)=squeeze(rbinned(:,respidx,attidx,:));
         nerfstd(:,:,1)=squeeze(rbinnedstd(:,respidx,attidx,:));
         
         for noiseidx=1:noisecount,
            [x,setidx]=sort(rand(size(rvalididx)));
            
            for eigidx=1:spacelimp,
               for bb=1:bincount,
                  sbbidx=find(tbinstim(setidx,eigidx)==bb);
                  if length(sbbidx)>0,
                     nerf(bb,eigidx,noiseidx+1)=...
                         mean(tr(sbbidx));
                     %mean(tr(sbbidx).*teigstim(sbbidx,eigidx))./...
                     %mean(teigstim(sbbidx,eigidx));
                     nerfstd(bb,eigidx,noiseidx+1)=...
                         std(tr(sbbidx)) ./ sqrt(length(sbbidx));
                  end
               end
            end
            if exist('BATQUEUEID') & BATQUEUEID>0,
               dbsetqueue(BATQUEUEID,respidx*1000+noiseidx);
            end
         end
         
         for eigidx=1:spacelimp,
            for noiseidx=1:noisecount+1,
               rm=nerf(:,eigidx,noiseidx);
               rs=nerfstd(:,eigidx,noiseidx);
               rs(find(rs==0))=1;
               m=rmean(respidx,attidx);
               fchi2(attidx,eigidx,noiseidx,respidx)=sum(abs(rm-m)./rs);
            end
            
            % measure p value of ERF modulation by comparing to
            % randomly generated ones
            nn=length(find(fchi2(attidx,eigidx,2:end,respidx)>=...
                           fchi2(attidx,eigidx,1,respidx)));               
            perf(attidx,eigidx,respidx)=(nn+1)./(noisecount+1);
         end
end
if 1
         rvalididx=find(~isnan(bresp(:,respidx,attidx)));
         nerf=zeros(bincount,spacelimp,noisecount+1);
         nerfstd=zeros(bincount,spacelimp,noisecount+1);
         tbinstim=binstim(rvalididx,:);
         teigstim=eigstim(rvalididx,:);
         tr=bresp(rvalididx,respidx,attidx);
         
         for eigidx=1:spacelimp,
            for bb=1:bincount,
               sbbidx=find(tbinstim(:,eigidx)==bb);
               if length(sbbidx)>0,
                  rbinned(bb,respidx,attidx,eigidx)=...
                      mean(tr(sbbidx));
                  %mean(tr(sbbidx).*teigstim(sbbidx,eigidx))./...
                  %mean(teigstim(sbbidx,eigidx));
                  rbinnedstd(bb,respidx,attidx,eigidx)=...
                      std(tr(sbbidx)) ./ sqrt(length(sbbidx));
               end
            end
         end
         
         nerf(:,:,1)=squeeze(rbinned(:,respidx,attidx,:));
         nerfstd(:,:,1)=squeeze(rbinnedstd(:,respidx,attidx,:));
         
         for noiseidx=1:noisecount,
            [x,setidx]=sort(rand(size(rvalididx)));
            
            for eigidx=1:spacelimp,
               for bb=1:bincount,
                  sbbidx=find(tbinstim(setidx,eigidx)==bb);
                  if length(sbbidx)>0,
                     nerf(bb,eigidx,noiseidx+1)=...
                         mean(tr(sbbidx));
                     %mean(tr(sbbidx).*teigstim(sbbidx,eigidx))./...
                     %mean(teigstim(sbbidx,eigidx));
                     nerfstd(bb,eigidx,noiseidx+1)=...
                         std(tr(sbbidx)) ./ sqrt(length(sbbidx));
                  end
               end
            end
            if exist('BATQUEUEID') & BATQUEUEID>0,
               dbsetqueue(BATQUEUEID,respidx*1000+noiseidx);
            end
         end
         
         for eigidx=1:spacelimp,
            for noiseidx=1:noisecount+1,
               rm=nerf(:,eigidx,noiseidx);
               rs=nerfstd(:,eigidx,noiseidx);
               rs(find(rs==0))=1;
               m=rmean(respidx,attidx);
               fchi2(attidx,eigidx,noiseidx,respidx)=sum(abs(rm-m)./rs);
            end
            
            % measure p value of ERF modulation by comparing to
            % randomly generated ones
            nn=length(find(fchi2(attidx,eigidx,2:end,respidx)>=...
                           fchi2(attidx,eigidx,1,respidx)));               
            perf(attidx,eigidx,respidx)=(nn+1)./(noisecount+1);
         end
end
if 0
         % compute random ERFs
         ntot=sum(rbincount(:,respidx,attidx,1));
         ntot=round(linspace(0,ntot,bincount+1));
         rset=cat(1,rbinset{:,respidx,attidx,1});
         for noiseidx=1:noisecount,
            [x,setidx]=sort(rand(size(rset)));
            
            % only test channels with >0 entries
            for bb=find(diff(ntot)>0),
               rtemp=rset(setidx((ntot(bb)+1):ntot(bb+1)));
               rm(bb)=mean(rtemp)-rmean(respidx,attidx);
               rs(bb)=std(rtemp) ./ sqrt(ntot(bb+1)-ntot(bb)+1);
            end
            
            m=0;
            rs(find(rs==0))=1;
               fchi2noise(attidx,noiseidx,respidx)=sum(abs(rm-m)./rs);
         end
         
         % compute ERFs for actual data
         %rbinned=zeros(bincount,respcount,attcount,spacelimp);
         %rbinnedstd=zeros(bincount,respcount,attcount,spacelimp);
         if ~isnan(rmean(respidx,attidx)),
            for eigidx=1:spacelimp,
               % compute mean and variance for each stimulus bin
               for bb=1:bincount,
                  if length(rbinset{bb,respidx,attidx,eigidx})>0,
                     rbinned(bb,respidx,attidx,eigidx)=...
                         mean(rbinset{bb,respidx,attidx,eigidx});
                     rbinnedstd(bb,respidx,attidx,eigidx)=...
                         std(rbinset{bb,respidx,attidx,eigidx}) ./ ...
                         sqrt(rbincount(bb,respidx,attidx,eigidx));
                  else
                     rbinned(bb,respidx,attidx,eigidx)=...
                         rmean(respidx,attidx);
                  end
               end
               
               rm=rbinned(:,respidx,attidx,eigidx);
               rs=rbinnedstd(:,respidx,attidx,eigidx);
               rs(find(rs==0))=1;
               m=rmean(respidx,attidx);
               fchi2act(attidx,eigidx,respidx)=sum(abs(rm-m)./rs);
               
               % measure p value of ERF modulation by comparing to
               % randomly generated ones
               nn=length(find(fchi2noise(attidx,:,respidx)>=...
                              fchi2act(attidx,eigidx,respidx)));               
               perf(attidx,eigidx,respidx)=(nn+1)./(noisecount+1);
            end
         end
end
      end

      if exist('BATQUEUEID') & BATQUEUEID>0,
         dbsetqueue(BATQUEUEID,respidx*1000);
      end
   end
   
   % 2.5 determine which spatial dimensions (PCs) show tuning
   
   % record original number of spatial dims tested
   startspacelimp=spacelimp;
   
   % find valid spatial dims
   % using erf finds more "valid" spatial bins.  is this better or worse?
   
   ptest=perf; % permute(pstrf,[3 1 2]) is worse?
   if min(ptest(:))>PFILT,
      PFILT=min(ptest(:))+PFILT;
   end
   spaceusep=find(sum(sum(ptest<PFILT,1),3)>0);
   spacelimp=length(spaceusep);
   
   fprintf('%d significant dims: ',spacelimp);
   for eigidx=spaceusep(:)',
      fprintf('%d ',eigidx);
   end
   fprintf(' (p<%.3f)\n',PFILT);
   
   % remove unwanted dims from existing variables
   %rbinset=reshape({rbinset{:,:,:,spaceusep}},bincount,respcount,...
   %                attcount,spacelimp);
   rbinned=rbinned(:,:,:,spaceusep);
   rbinnedstd=rbinnedstd(:,:,:,spaceusep);
   rbincount=rbincount(:,:,:,spaceusep);
   sbin=sbin(:,spaceusep);
   mpbinned=mpbinned(:,spaceusep);
   eigfpatches=eigfpatches(:,spaceusep);
   eigmpatches=eigmpatches(:,spaceusep);
   eigstim=eigstim(:,spaceusep);
   binstim=binstim(:,spaceusep);
   
   % display ERF for each PC where it shows significant modulation
   figure(2);
   clf
   respidx=min([respcount 7]);
   smark={'k','g','c','b','r'};
   %rm=reshape(z.brmean(:,:,2:end),respcount,attcount-1,noisesets);
   %rmall=squeeze(mean(rm(respidx,:,1)));
   rmall=0;
   rmin=min(min(min(min(rbinned(:,respidx,2:attcount,:)))))-rmall;
   rmax=max(max(max(max(rbinned(:,respidx,2:attcount,:)))))-rmall;
   if (rmax-rmin)==0,
      rmax=rmin+1;
   end
   spaceshow=spacelimp;
   %spaceshow=15;
   rowcount=ceil(spaceshow/5);
   for eigidx=1:spaceshow,
      subplot(rowcount,ceil(spaceshow/rowcount),eigidx);
      
      if 0,
         imagesc(squeeze(rbinned(:,respidx,1:attcount,eigidx))',[rmin,rmax]);
         hold on
         for attidx=1:attcount,
            plot(mpbinned(attidx,eigidx),attidx,'kx');
            plot(mpbinned(attidx,eigidx),attidx,'wo');
         end
         hold off
         
         axis off;
      else
         hold on
         labels={'A','1','2','3','4'};
         ll=zeros(attcount,1);
         attshowset=1:attcount;
         for attidx=attshowset,
            r=rbinned(:,respidx,attidx,eigidx);
            %r=rbinned(:,respidx,attidx,eigidx)-...
            %  mean(rbinned(:,respidx,attidx,eigidx));
            a=errorbar(sbin(:,eigidx),r,...
                       squeeze(rbinnedstd(:,respidx,attidx,eigidx)).*2,...
                       [smark{attidx},'-']);
            %set(a,'LineWidth',2);
            ll(attidx)=a(1);
            plot(sbin(mpbinned(attidx,eigidx),eigidx),...
                 r(mpbinned(attidx,eigidx)),'kx');
         end
         hold off
         axis([sbin(1,eigidx),sbin(end,eigidx),rmin,rmax]);
         axis([sbin(1,eigidx),sbin(end,eigidx),rmin,rmax]);
         if eigidx==spaceshow,
            legend(ll(find(ll)),labels{find(ll)});
         end
      end
      title(sprintf('%s: pc%d',cellid,eigidx));
   end
   colormap(hot);
   set(gcf,'PaperOrientation','Landscape','PaperPosition',[0.25 0.25 10.5 8]);
   %set(gcf,'PaperOrientation','Landscape','PaperPosition',[0.25 2 10.5 4]);
   
   if 0 & strcmp(cellid,'m0076'),
      keyboard;
   end
   
   if 0,
      ep=[1 2; 1 3; 1 4; 2 3; 2 4; 3 4]+30;
      epcount=size(ep,1);
      test2d=zeros(bincount,bincount,attcount,epcount);
      test2dse=zeros(bincount,bincount,attcount,epcount);
      for epidx=1:epcount,
         e1=ep(epidx,1);
         e2=ep(epidx,2);
         for attidx=1:attcount,
            for b1=1:bincount,
               for b2=1:bincount,
                  tidx=find(binstim(:,e1)==b1 & binstim(:,e2)==b2 & ...
                            ~isnan(bresp(:,1,attidx)));
                  if length(tidx>0),
                     test2d(b1,b2,attidx,epidx)=mean(bresp(tidx,1,attidx));
                     test2dse(b1,b2,attidx,epidx)=...
                         std(bresp(tidx,1,attidx))./length(tidx);
                  end
               end
            end
         end
      end
      
      figure(4);
      clf
      for epidx=1:epcount,
         e1=ep(epidx,1);
         e2=ep(epidx,2);
         rmin=min(min(min(test2d(:,:,:,epidx))));
         rmax=max(max(max(test2d(:,:,:,epidx))));
         for attidx=1:attcount,
            subplot(epcount,attcount,(epidx-1)*attcount+attidx);
            %[xx,yy]=meshgrid(sbin(:,e1),sbin(:,e2));
            %surf(xx,yy,test2d);
            imagesc(sbin(:,e1),sbin(:,e2),test2d(:,:,attidx,epidx),[rmin,rmax]);
            hold on
            plot(sbin(mpbinned(attidx,e1),e1),sbin(mpbinned(attidx,e2),e2),'kx');
            plot(sbin(mpbinned(attidx,e1),e1),sbin(mpbinned(attidx,e2),e2),'wo');
            hold off
            axis square
            xlabel(sprintf('pc%d',e1));
            ylabel(sprintf('pc%d',e2));
         end
         colormap(hot);
      end
   end
   
   % 3. determine significantly attention modulated spatial channels
   % (independent of attention)
   disp('3. Calculating ERF attention-modulation T2 and significance...');
   %aT2act=zeros(attdocount,spacelimp,respcount,2);
   %aT2noise=zeros(attdocount,noisecount,spacelimp,respcount,2);
   aT2=zeros(attdocount,spacelimp,noisecount+1,respcount,3);
   perfatt0=ones(attdocount,spacelimp,respcount);
   perfatt=ones(attdocount,spacelimp,respcount);
   perfatt1=ones(attdocount,spacelimp,respcount);
   
   for respidx=1:respcount,
      fprintf('respidx=%d: paircomp ',respidx);
      
      for pidx=1:attdocount,
         a1=attpairs(pidx,1)+1;
         a2=attpairs(pidx,2)+1;
         fprintf('%d v %d  ',a1,a2);
         
if 1         
         r1valididx=find(~isnan(bresp(:,respidx,a1)));
         r2valididx=find(~isnan(bresp(:,respidx,a2)));
         rvalididx=[r1valididx; r2valididx];
         nerf=zeros(bincount,spacelimp,noisecount+1,2);
         nerfstd=zeros(bincount,spacelimp,noisecount+1,2);
         nmean=zeros(noisecount+1,2);
         
         tbinstim=binstim(rvalididx,:);
         teigstim=eigstim(rvalididx,:);
         tr=[bresp(r1valididx,respidx,a1);
             bresp(r2valididx,respidx,a2)];
         
         nerf(:,:,1,1)=squeeze(rbinned(:,respidx,a1,:));
         nerfstd(:,:,1,1)=squeeze(rbinnedstd(:,respidx,a1,:));
         nerf(:,:,1,2)=squeeze(rbinned(:,respidx,a2,:));
         nerfstd(:,:,1,2)=squeeze(rbinnedstd(:,respidx,a2,:));
         nmean(1,:)=rmean(respidx,[a1 a2]);
         
         for noiseidx=1:noisecount,
            [x,setidx]=sort(rand(size(rvalididx)));
            n1idx=setidx(1:length(r1valididx));
            n2idx=setidx(length(r1valididx)+1:end);
            
            nmean(noiseidx+1,1)=mean(tr(n1idx));
            nmean(noiseidx+1,2)=mean(tr(n2idx));
            
            for eigidx=1:spacelimp,
               for bb=1:bincount,
                  sbbidx=find(tbinstim(n1idx,eigidx)==bb);
                  if length(sbbidx)>0,
                     nerf(bb,eigidx,noiseidx+1,1)=...
                         mean(tr(n1idx(sbbidx)));
                        %mean(tr(sbbidx).*teigstim(sbbidx,eigidx))./...
                        %mean(teigstim(sbbidx,eigidx));
                     nerfstd(bb,eigidx,noiseidx+1,1)=...
                         std(tr(n1idx(sbbidx)))./sqrt(length(sbbidx));
                  end
                  
                  sbbidx=find(tbinstim(n2idx,eigidx)==bb);
                  if length(sbbidx)>0,
                     nerf(bb,eigidx,noiseidx+1,2)=...
                         mean(tr(n2idx(sbbidx)));
                     %mean(tr(sbbidx).*teigstim(sbbidx,eigidx))./...
                     %mean(teigstim(sbbidx,eigidx));
                     nerfstd(bb,eigidx,noiseidx+1,2)=...
                         std(tr(n2idx(sbbidx)))./sqrt(length(sbbidx));
                  end
               end
            end
            if exist('BATQUEUEID') & BATQUEUEID>0,
               dbsetqueue(BATQUEUEID,respidx*1000+pidx*100+noiseidx);
            end
         end
         
         %sad=diag(sa);
         sad=ones(128,1);
         for eigidx=1:spacelimp,
            for noiseidx=1:noisecount+1,
               rm=diff(nerf(:,1:eigidx,noiseidx,:),1,4) .* ...
                  sqrt(repmat(sad(1:eigidx)',[bincount,1]));
               rs=mean(nerfstd(:,1:eigidx,noiseidx,:),4);
               rs(find(rs==0))=1;
               
               aT2(pidx,eigidx,noiseidx,respidx,1)=sum(abs(rm(:)./rs(:)));
               
               m=diff(nmean(noiseidx,:)) .* ...
                  sqrt(repmat(sad(1:eigidx)',[bincount,1]));
               aT2(pidx,eigidx,noiseidx,respidx,2)=...
                   sum(abs((rm(:)-m(:))./rs(:)));
               
               rml=squeeze(rm(end,:));
               ml=squeeze(m(end,:));
               rsl=squeeze(rs(end,:));
               aT2(pidx,eigidx,noiseidx,respidx,3)=...
                   sum(abs((rml(:)-ml(:))./rsl(:)));
            end
            
            % measure p value of ERF modulation by comparing to
            % randomly generated ones
            nn=length(find(aT2(pidx,eigidx,2:end,respidx,1)>=...
                           aT2(pidx,eigidx,1,respidx,1)));               
            perfatt0(pidx,eigidx,respidx)=(nn+1)./(noisecount+1);
            
            % measure p value of ERF modulation by comparing to
            % randomly generated ones
            nn=length(find(aT2(pidx,eigidx,2:end,respidx,2)>=...
                           aT2(pidx,eigidx,1,respidx,2)));               
            perfatt(pidx,eigidx,respidx)=(nn+1)./(noisecount+1);
            
            % measure p value of ERF for individual eigs to see
            % where it's maximal
            nn=length(find(aT2(pidx,eigidx,2:end,respidx,3)>=...
                           aT2(pidx,eigidx,1,respidx,3)));               
            perfatt1(pidx,eigidx,respidx)=(nn+1)./(noisecount+1);
            
         end
end
         
if 0
         % T2 for actual data
         for eigidx=1:spacelimp,
            
            % actual attentional conditions
            x1=squeeze(rbinned(:,respidx,a1,1:eigidx));
            x2=squeeze(rbinned(:,respidx,a2,1:eigidx));
            x1std=squeeze(rbinnedstd(:,respidx,a1,1:eigidx));
            x2std=squeeze(rbinnedstd(:,respidx,a2,1:eigidx));
            
            xstd=(x1std+x2std)./2;
            xstd(find(xstd==0))=1;
            aT2act(pidx,eigidx,respidx,1)=...
                sum(abs(x1(:)-x2(:))./xstd(:),1);
            
            y1=x1(:,1:eigidx)-rmean(respidx,a1);
            y2=x2(:,1:eigidx)-rmean(respidx,a2);
            ystd=xstd;
            aT2act(pidx,eigidx,respidx,2)=...
                sum(abs(y1(:)-y2(:))./ystd(:),1);
         end
         
         % T2 for noise (shuffled) conditions
         for noiseidx=1:noisecount,
            x1=zeros(bincount,spacelimp);
            x2=zeros(bincount,spacelimp);
            x1std=zeros(bincount,spacelimp);
            x2std=zeros(bincount,spacelimp);
            xstd=zeros(bincount,spacelimp);
            m1=zeros(1,spacelimp);
            m2=zeros(1,spacelimp);
            n1=zeros(1,spacelimp);
            n2=zeros(1,spacelimp);
            for eigidx=1:spacelimp,
               for bb=1:bincount,
                  ntot=squeeze(rbincount(bb,respidx,[a1 a2],eigidx));
                  rset=cat(1,rbinset{bb,respidx,[a1 a2],eigidx});
                  
                  [x,setidx]=sort(rand(size(rset)));
                  
                  if ntot(1)>0,
                     r1=rset(setidx(1:ntot(1)));
                     m1(eigidx)=m1(eigidx)+sum(r1);
                     n1(eigidx)=n1(eigidx)+length(r1);
                     x1(bb,eigidx)=mean(r1);
                     x1std(bb,eigidx)=std(r1) ./ sqrt(ntot(1));
                  end
                  if ntot(2)>0,
                     r2=rset(setidx((ntot(1)+1):end));
                     m2(eigidx)=m2(eigidx)+sum(r2);
                     n2(eigidx)=n2(eigidx)+length(r2);
                     x2(bb,eigidx)=mean(r2);
                     x2std(bb,eigidx)=std(r2) ./ sqrt(ntot(2));
                  end
               end
               
               % avg std across both att states fo T2
               xstd(:,eigidx)=(x1std(:,eigidx)+x2std(:,eigidx))./2;
               xstd(find(xstd(:,eigidx)==0),eigidx)=1;
               
               % mean response for each attstate/PC
               m1(eigidx)=m1(eigidx)./n1(eigidx);
               m2(eigidx)=m2(eigidx)./n2(eigidx);
               
               aT2noise(pidx,noiseidx,eigidx,respidx,1)=...
                   sum(sum(abs(x1(:,1:eigidx)-x2(:,1:eigidx))./...
                       xstd(:,1:eigidx),1));
               %./ ...
               %  (sqrt(sum(sum(x1(:,1:eigidx).^2./xstd(:,1:eigidx),1))).* ...
               %     sqrt(sum(sum(x2(:,1:eigidx).^2./xstd(:,1:eigidx),1))));
               
               y1=x1(:,1:eigidx)-repmat(m1(1:eigidx),[bincount,1]);
               y2=x2(:,1:eigidx)-repmat(m2(1:eigidx),[bincount,1]);
               ystd=xstd(:,1:eigidx);
               
               aT2noise(pidx,noiseidx,eigidx,respidx,2)=...
                   sum(abs(y1(:)-y2(:))./ystd(:));
               %./ ...
               %   (sqrt(sum(y1(:).^2./ystd(:))).* ...
               %   sqrt(sum(y2(:).^2./ystd(:))));
            end
            if exist('BATQUEUEID') & BATQUEUEID>0,
               dbsetqueue(BATQUEUEID,300+respidx*10+pidx);
            end
         end
         
         % now that all the noise T2's have been generated, compute
         % p for out to each PC
         for eigidx=1:spacelimp,
            % compute significance
            nn=length(find(aT2noise(pidx,:,eigidx,respidx,1)>...
                           aT2act(pidx,eigidx,respidx,1)));
            perfatt0(pidx,eigidx,respidx)=(nn+1)./(noisecount+1);
            
            % compute significance
            nn=length(find(aT2noise(pidx,:,eigidx,respidx,2)>=...
                           aT2act(pidx,eigidx,respidx,2)));
            perfatt(pidx,eigidx,respidx)=(nn+1)./(noisecount+1);
         end
end
         if exist('BATQUEUEID') & BATQUEUEID>0,
            dbsetqueue(BATQUEUEID,300+respidx*10+pidx);
         end
      end
      fprintf('\n');
   end
   
   %keyboard
   
   % assess the results. for each pair, what dims were
   % significantly different?
   %
   % pstrf, perf - attention indep modulation
   %  (attcount x spacelim0 x respcount) 
   % pstrfatt - attention modulation of strf
   % perfatt0 - attention modulation of erf at all (with mean)
   % perfatt - attention modulation of erf tuning (mean subtracted)
   %  (attdocount x spacelim x respcount)
   %
   figure(3);
   clf
   for eigidx=1:spacelimp,
      subplot(ceil(spacelimp/6),6,eigidx);
      pdata=[squeeze(perfatt0(:,eigidx,:)); ...
             squeeze(perfatt(:,eigidx,:))]';
      imagesc(pdata,[0 1]);
      
      hold on
      for respidx=1:respcount,
         tt=find(pdata(respidx,:)<PATT);
         if length(tt)>0,
            plot(tt,respidx,'bx');
         end
      end
      hold off
      
      if eigidx==1,
         title(sprintf('%s eigidx=%d',cellid,eigidx));
      else
         title(sprintf('eigidx=%d',eigidx));
      end
   end
   colormap(flipud(hot))
   set(gcf,'PaperOrientation','Landscape','PaperPosition',[0.25 0.25 10.5 8]);
   
   
   if 0,
   torder=get(0,'DefaultAxesColorOrder');
   set(0,'DefaultAxesColorOrder',[0 0 0],...
         'DefaultAxesLineStyleOrder','-|--|:|-.');
   figure(4);
   clf
   subplot(2,1,1);
   % find tuned dims
   if min([pstrf(:); perf(:)])>PFILT,
      PFILT=min(pstrf(:));
   end
   pdata=[sum(sum(pstrf<PFILT,2),3) sum(sum(perf<PFILT,1),3)'];
   plot(pdata);
   legend('strf','erf');
   title(sprintf('%s stimulus tuning',cellid));
   ylabel(sprintf('n(att/resp) p<%.3f',PFILT));
   xlabel('PC idx');
   
   subplot(2,1,2);
   % find att mod dims
   pdata=[sum(sum(pstrfatt<PATT,1),3)' ...
          sum(sum(perfatt0<PATT,1),3)' ...
          sum(sum(perfatt<PATT,1),3)'];
   plot(spaceusep,pdata);
   legend('strf','erf0','erf');
   title(sprintf('%s attention modulation',cellid));
   ylabel(sprintf('n(att/resp) mod p<%.3f',PATT));
   xlabel('PC idx');
   set(0,'DefaultAxesColorOrder',torder);
   
   disp('just finished kerncomp3');
   drawnow
   end
end

targsim=zeros(attdocount,spacelimp,2);
attcc=zeros(spacelimp,attcount,3);
p=0;

% clear bresp
clear bstim x1 x2 r2 r2std
clear data rbinset y1 y2 ystd *idx


return

