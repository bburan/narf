% kerncomp.m : compare kernels from different attentional states to
% see how similar they are

disp('kerncomp4.m : back to a more traditional strf');

% parms for XC / decorr
if ~exist('maxlag','var'),
   %maxlag=[-10 15];
   maxlag=[0 0];
end
if ~exist('meansub','var'),
   meansub=1;
end
if ~exist('boundary','var'),
   boundary='zero';
end
if ~exist('mineigs','var'),
   mineigs=[];
end
% movstep: size of movie segments to send to movxc -- want this to
% be as large as possible without straining computer's memory
if ~exist('movstep','var'),
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

spacelimp=20;
PFILT=0.010;  % p value (any attention) to be considered stimulus-modulated
PATT=0.008;   % p<0.008 any pair has p<ATTA attention modulation
PATTA=0.05;   % p for attention modulation

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
   
   % optionally only do analysis on att states with biggest mean
   % difference in response
   if 0,
      respidx=min([size(bresp,2) 7]);
      nresp=squeeze(nansum(~isnan(bresp(:,respidx,:))));
      mresp=squeeze(nanmean(bresp(:,respidx,:)));
      attmin=find(mresp(2:end)==min(mresp(2:end)))+1;
      attmax=find(mresp(2:end)==max(mresp(2:end)))+1;
      bresp=bresp(:,:,[1 attmin attmax]);
      %attlong=find(nresp>100);
      %bresp=bresp(:,:,attlong);
   end
   
   % figure out size and mean of response
   rsize=size(bresp);
   blen=rsize(1);
   respcount=rsize(2);
   attcount=rsize(3);
   
   % adjust attcount to actual number of attention states
   if fidx==1,
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
   
   DAMPFACTOR=0.8;
   DODAMPING=0;
   bootcount=20;
   spacelims=128;
   %spacelims=35;
   xccount=2;
   spaceuses=1:spacelims;
   predspace=spaceuses;   % if active cutoff selection is on, this
                          % will be replaced after the first iteration
   
   % global gain for each attentional kernel
   a1=zeros(respcount,attcount,noisecount+1,bootcount,xccount);
   % dc offset for each attentional state
   c1=zeros(respcount,attcount,noisecount+1,bootcount,xccount);
   % full (locally modulated) kernel for each attentional state
   SR=zeros(spacecount,respcount,attcount,noisecount+1,bootcount);
   Hpc=zeros(spacelims,respcount,attcount,noisecount+1,bootcount);
   
   % placeholders
   lSR=zeros(spacecount,respcount,attcount,1,bootcount);
   lHpc=zeros(spacelims,respcount,attcount,1,bootcount);
   
   % prediction scores for real and randomized attention
   % conditions, exploratory and confirmatory response data
   paircount=(attcount-1)*(attcount-2)/2;
   xc=zeros(noisecount+1,4,xccount);
   xcp=zeros(noisecount+1,4,xccount,paircount);
   
   % generate PC domain basis set for projecting all H's
   aokidx=find(sum(~isnan(bresp(:,end,2:end)),3));
   stim=bstim(aokidx,:);
   sm=mean(stim,1);
   stim=stim-repmat(sm,[size(stim,1),1]);
   sSAfull=stim'*stim ./ length(aokidx);
   [ua,sa,va]=svd(sSAfull);
   eigstim=bstim*ua;
   %keyboard
   
   % copy resp and subtract mean from each attentional state in
   % order to support randomization over gain/local effects only
   % (ie, randomized kernels have appropriate dc shifts)
   brespnodc=bresp;
   for attidx=2:attcount,
      for respidx=1:respcount,
         aokidx=find(~isnan(brespnodc(:,respidx,attidx)));
         brespnodc(aokidx,respidx,attidx)=...
             brespnodc(aokidx,respidx,attidx)-...
             mean(brespnodc(aokidx,respidx,attidx));
         brespnodc(aokidx,respidx,1)=brespnodc(aokidx,respidx,attidx);
      end
   end
   meanadjust=bresp(:,:,1)-brespnodc(:,:,1);
   
   % figure out number of valid fixations in each attentional condition
   ncount=squeeze(sum(~isnan(bresp(:,end,:)),1));
   cumncount=[0;cumsum(ncount(2:end))];
   nrange=find(sum(~isnan(bresp(:,end,2:end)),3));
   
   % choose a random set of fixations for each bootstrap. this should
   % avoid bias in the non-randomized noiseidx=1
   vidx=nrange;
   vcount=length(vidx);
   [vv,ss]=sort(rand(size(vidx)));
   vidx=vidx(ss);
   
   fprintf('DODAMPING=%d  DAMPFACTOR=%.2f  length(spaceuses)=%d\n',...
           DODAMPING,DAMPFACTOR,length(spaceuses));
   
   % loop 2 times, once on raw resp, once on resp with att-dc
   % already taken out.
   for xcidx=1:xccount,
      for noiseidx=1:noisecount+1,
         fprintf('noiseidx=%d',noiseidx);
         
         resp=ones(size(bresp))*nan;
         if xcidx==1,
            resp(nrange,:,1)=bresp(nrange,:,1);
         else
            resp(nrange,:,1)=brespnodc(nrange,:,1);
         end
         if noiseidx==1,
            if xcidx==1,
               resp(:,:,2:end)=bresp(:,:,2:end);
            else
               resp(:,:,2:end)=brespnodc(:,:,2:end);
            end
         else
            % pick random atts
            [sn,nidx]=sort(rand(size(nrange)));
            for attidx=2:attcount,
               nn=nrange(nidx(cumncount(attidx-1)+1:cumncount(attidx)));
               if xcidx==1,
                  resp(nn,:,attidx)=bresp(nn,:,1);
               else
                  resp(nn,:,attidx)=brespnodc(nn,:,1);
               end
            end
         end
         %keyboard
         smean=zeros(spacecount,respcount,attcount,bootcount);
         
         for bootidx=1:bootcount;
            fprintf('.');
            
            % generate exploratory response matrix that is resp
            % with cnf data removed
            eresp=resp;
            ppidx=vidx(round((bootidx-1)/bootcount*vcount+1):...
                       round(bootidx/bootcount*vcount));
            eresp(ppidx,:,:)=nan;
            ermean=reshape(nanmean(eresp),respcount,attcount);
            ermean(find(isnan(ermean)))=0;
            
            % find DC offset for this bootstrap
            c1(:,:,noiseidx,bootidx,xcidx)=ermean;
            
            for attidx=1:attcount,
               
               % find spatial correlation for stim in attidx
               
               % non-zero frames for this att/resample
               aokidx=find(~isnan(eresp(:,end,attidx)));
               
               % set up B to transform strf from attention state idx
               % into PC domain of all att, decorrelating as needed.
               if length(aokidx)>0,
                  
                  % ssa=u*s*v'
                  % ssa_pc=ua'* (u * s * v') * ua
                  % ssa^-1=v*s^-1*u'
                  % ssa_pc^-1 = ua'* v* s^-1 * u'* ua
                  if attidx==1,
                     % basically just have global correlations.
                     sst=[1./sa(find(sa>0)); ...
                          zeros(length(find(diag(sa)==0)),1)];
                     B=diag(sst(spaceuses)) * ua(:,spaceuses)';
                  else
                     stim=bstim(aokidx,:);
                     sm=mean(stim,1);
                     stim=stim-repmat(sm,[size(stim,1),1]);
                     tSA=stim'*stim ./ length(aokidx);
                     
                     if 0,
                        [ut,st,vt]=svd(tSA);
                        sst=[1./st(find(st>0)); ...
                             zeros(length(find(diag(st)==0)),1)];
                        B=ua(:,spaceuses)' * vt * diag(sst) * ut';
                     else
                        % alt, don't decorr, better SNR?????
                        st=diag(ua' * tSA *ua);
                        sst=[1./st(find(st>0)); ...
                             zeros(length(find(st<=0)),1)];
                        B=diag(sst(spaceuses)) * ua(:,spaceuses)';
                     end
                  end
                  
                  for respidx=1:respcount,
                     okidx=find(~isnan(eresp(:,respidx,attidx)));
                     slen=length(okidx);
                     if slen>0,
                        r=eresp(okidx,respidx,attidx)-ermean(respidx,attidx);
                        stim=bstim(okidx,:);
                        smean(:,respidx,attidx,bootidx)=mean(stim,1)';
                        if attidx==1 & noiseidx>1,
                           % attidx==1 is always the same for each noiseidx
                           lSR(:,respidx,attidx,1,bootidx)=...
                               SR(:,respidx,attidx,1,bootidx);
                           lHpc(:,respidx,attidx,1,bootidx)=...
                               Hpc(:,respidx,attidx,1,bootidx);
                        else
                           stim=stim-...
                                repmat(smean(:,respidx,attidx,bootidx)',...
                                       [slen,1]);
                           
                           lSR(:,respidx,attidx,1,bootidx)=stim'*r./slen;
                           lHpc(:,respidx,attidx,1,bootidx)=...
                               B*lSR(:,respidx,attidx,1,bootidx);
                        end
                        
                        % figure out optimal global gain for att all
                        % to pred current attidx responses
                        %estim=stim*ua(:,spaceuses);
                        estim=eigstim(okidx,spaceuses)-...
                              repmat(smean(:,respidx,attidx,bootidx)' ...
                                     * ua(:,spaceuses), [slen 1]);
                        r1=estim * lHpc(spaceuses,respidx,1,1,bootidx);
                        d1=sum(r1.^2);
                        if d1>0,
                           a1(respidx,attidx,noiseidx,bootidx,xcidx)=...
                               sum(r.*r1)./d1;
                        else
                           fprintf('z');
                           a1(respidx,attidx,noiseidx,bootidx,xcidx)=1;
                        end
                     end
                  end
               end
               if exist('BATQUEUEID','var') & BATQUEUEID>0,
                  dbsetqueue(BATQUEUEID,noiseidx*10+attidx);
               end
            end
         end
         
         % do prediction test noatt/att comparison
         respidx=1;
         ract=bresp(:,respidx,1);
         
         if xcidx==1 & noiseidx==1,
            r0=zeros(blen,1);
            
            Hm=mean(lHpc(:,respidx,1,1,:),5);
            Hs=std(lHpc(:,respidx,1,1,:),1,5) .* sqrt(bootcount);
            Hs(find(Hs==0))=1;
            
            ddr=0.3:0.1:1.5;
            xct=zeros(length(ddr),1);
            agoodidx=find(sum(~isnan(bresp(:,end,2:attcount)),3));
            
            for dd=1:length(ddr),
               Hdamp=(abs(Hm)./(Hs.*ddr(dd)));
               Hdamp=(1-Hdamp.^(-2));
               Hdamp=Hdamp.*(Hdamp>0);
               Hdamp(find(isnan(Hdamp)))=0;
               
               for bootidx=1:bootcount,
                  a0idx=vidx(round((bootidx-1)/bootcount*vcount+1):...
                             round(bootidx/bootcount*vcount));
                  
                  % pull out stimulus
                  estim=eigstim(a0idx,spaceuses)- ...
                     repmat(smean(:,respidx,1,bootidx)' * ...
                            ua(:,spaceuses), [length(a0idx) 1]);
                  
                  % prepare kernels for prediction with damping
                  H0=lHpc(predspace,respidx,1,1,bootidx).*Hdamp(predspace,1);
                  
                  % no attention prediction
                  r0(a0idx)=estim * H0;
                  r0(a0idx)=r0(a0idx) + c1(respidx,1,noiseidx,bootidx,xcidx);
               end
               xct(dd)=xcov(r0(agoodidx),ract(agoodidx),0,'coeff');
            end
            
            % generate optimal Hdamp for the remainder of predictions
            ddmax=find(xct==max(xct));
            Hdamp=(abs(Hm)./(Hs.*ddr(ddmax)));
            Hdamp=(1-Hdamp.^(-2));
            Hdamp=Hdamp.*(Hdamp>0);
            Hdamp(find(isnan(Hdamp)))=0;
            Hdamp=repmat(Hdamp,[1 attcount]);
            
            fprintf('\nmaximum DAMPFACTOR=%.1f\n',ddr(ddmax));
         end

         % determine optimal number of eigenvectors to include in
         % kernel using all-att condition
         %keyboard
         if ~DODAMPING & xcidx==1 & noiseidx==1,
            r0=zeros(blen,spacelims);
            for bootidx=1:bootcount,
               a0idx=vidx(round((bootidx-1)/bootcount*vcount+1):...
                       round(bootidx/bootcount*vcount));
            
               % pull out stimulus
               estim=eigstim(a0idx,spaceuses)- ...
                     repmat(smean(:,respidx,1,bootidx)' * ...
                            ua(:,spaceuses), [length(a0idx) 1]);
            
               % prepare kernels for prediction
               % no damping
               H0=lHpc(spaceuses,respidx,1,1,bootidx);
               
               % no attention prediction
               r0(a0idx,:)=estim .* repmat(H0',[length(a0idx),1]);
               r0(a0idx)=r0(a0idx) + c1(respidx,1,noiseidx,bootidx,xcidx);
            end
            
            agoodidx=find(sum(~isnan(bresp(:,end,2:attcount)),3));
            xct=zeros(spacelims,1);
            for ss=1:spacelims,
               xct(ss)=xcov(sum(r0(agoodidx,1:ss),2),ract(agoodidx),...
                                      0,'coeff');
            end
            
            % smooth with boxcar to reduce local maxima
            xct=conv2(xct,[1 1 1]','same');
            emax=find(xct==max(xct));
            predlims=emax;
            predspace=1:predlims;
            
            fprintf('maximum pred at predlims=%d eigs\n', ...
                    predlims);
         end
         
         % set up vectors to hold preds from different models
         r0=zeros(blen,1);
         rnoatt=zeros(blen,1);
         rdcatt=zeros(blen,1);
         rglobalatt=zeros(blen,1);
         rfullatt=zeros(blen,1);
         
         % compute contribution of each bootstrap component to the preds.
         for bootidx=1:bootcount,
            if DODAMPING==2,
               Hm=mean(lHpc(:,respidx,:,1,...
                            [1:bootidx-1 bootidx+1:bootcount]),5);
               Hs=std(lHpc(:,respidx,:,1,...
                           [1:bootidx-1 bootidx+1:bootcount]),1,5)...
                  .*sqrt(bootcount-1);
               
               % slightly cheating method for generating Hdamp.
               %Hm=mean(lHpc(:,:,:,1,:),5);
               %Hs=std(lHpc(:,:,:,1,:),1,5).*sqrt(bootcount);
               %Hs=repmat(std(lHpc(:,:,1,1,:),1,5).*sqrt(bootcount),
               %     [1 1 attcount]);
               
               %
               % shrinkage filter. is this better? seems to help a lot!
               % 
               Hs(find(Hs==0))=1;
               Hdamp=(abs(Hm)./(Hs.*DAMPFACTOR));
               Hdamp=(1-Hdamp.^(-2));
               Hdamp=Hdamp.*(Hdamp>0);
               Hdamp(find(isnan(Hdamp)))=0;
               
               % rescale to keep power of pred response about right.
               %for attidx=1:attcount,
               %   if sum(Hdamp(:,attidx))>0,
               %      Hdamp(:,attidx)=Hdamp(:,attidx)./ ...
               %          (sum(Hm(:,attidx).^2.*Hdamp(:,attidx).*diag(sa)) ...
               %           ./ sum(Hm(:,attidx).^2.*diag(sa)));
               %   end
               %end
            end
            
            % moved bootidx loop earlier so that Hdamp is not
            % biased by confirmatory set.  is this all i need to do
            % to make it legitimately unbiased?
            
            % find cnf fixations for this bootstrap
            a0idx=vidx(round((bootidx-1)/bootcount*vcount+1):...
                       round(bootidx/bootcount*vcount));
            
            % pull out stimulus
            %stim=bstim(a0idx,:);
            %stim=stim-repmat(smean(:,respidx,1,bootidx)',...
            %                 [size(stim,1),1]);
            %estim=stim*ua(:,predspace);
            estim=eigstim(a0idx,predspace)- ...
                  repmat(smean(:,respidx,1,bootidx)' * ...
                         ua(:,predspace), [length(a0idx) 1]);
            
            % prepare kernels for prediction
            if ~DODAMPING,
               % no damping
               H0=lHpc(predspace,respidx,1,1,bootidx);
               H=squeeze(lHpc(predspace,respidx,2:end,1,bootidx));
            else
               % damping on
               H0=lHpc(predspace,respidx,1,1,bootidx).*Hdamp(predspace,1);
               H=squeeze(lHpc(predspace,respidx,2:end,1,bootidx)) .* ...
                 Hdamp(predspace,2:end);
            end
            
            % no attention prediction
            r0(a0idx)=estim * H0;
            rnoatt(a0idx)=r0(a0idx) + c1(respidx,1,noiseidx,bootidx,xcidx);
            
            for attidx=2:attcount,
               % find prediction fixations for this bootstrap and attidx
               aidx=a0idx(find(~isnan(resp(a0idx,respidx,attidx))));
               
               % dc offset prediction
               rdcatt(aidx)=r0(aidx) + ...
                   c1(respidx,attidx,noiseidx,bootidx,xcidx);
               
               % global gain prediction
               rglobalatt(aidx)=r0(aidx) .* ...
                   a1(respidx,attidx,noiseidx,bootidx,xcidx) + ...
                   c1(respidx,attidx,noiseidx,bootidx,xcidx);
               
               % local gain prediction
               estim=eigstim(aidx,predspace)-...
                     repmat(smean(:,respidx,attidx,bootidx)' ...
                            * ua(:,predspace), [length(aidx) 1]);
               rfullatt(aidx)=estim * H(:,attidx-1) + ...
                   c1(respidx,attidx,noiseidx,bootidx,xcidx);
            end
         end
         
         % any attention state
         agoodidx=find(sum(~isnan(bresp(:,end,2:attcount)),3));
         
         % count how many PCs are significantly non-zero for each
         % attentional state (for actual local attention strfs)
         if noiseidx==1 & xcidx==1,
            nondampedcount=squeeze(sum(Hdamp(:,respidx,:)>0,1));
         end
         
         % for local attention only choose local atttention states
         % that had any significant power in their STRF
         if DODAMPING==2,
            goodatt=find(nondampedcount(2:end)>0)'+1;
         else
            goodatt=2:attcount; % old way, no excluded attidx's
         end
         alocalidx=find(sum(~isnan(resp(:,end,goodatt)),3));
         
         if xcidx==2,
            % add back dc shifts to make xc scores comparable
            rnoatt=rnoatt+meanadjust;
            rdcatt=rdcatt+meanadjust;
            rglobalatt=rglobalatt+meanadjust;
            rfullatt=rfullatt+meanadjust;
         end
         
         % rectify?  doesn't seem to help worth dog poo poo
         if 0,
            rnoatt(find(rnoatt<0))=0;
            rdcatt(find(rdcatt<0))=0;
            rglobalatt(find(rglobalatt<0))=0;
            rfullatt(find(rfullatt<0))=0;
         end
         
         if length(agoodidx)>0,
            xc(noiseidx,1,xcidx)=xcov(rnoatt(agoodidx),ract(agoodidx),...
                                      0,'coeff');
            xc(noiseidx,2,xcidx)=xcov(rdcatt(agoodidx),ract(agoodidx),...
                                      0,'coeff');
            xc(noiseidx,3,xcidx)=xcov(rglobalatt(agoodidx),...
                                      ract(agoodidx),0,'coeff');
            if length(alocalidx)>0,
               xc(noiseidx,4,xcidx)=xcov(rfullatt(alocalidx),...
                                         ract(alocalidx),0,'coeff');
            end
            fprintf('xc%d: %.3f %.3f %.3f %.3f\n',xcidx,xc(noiseidx,:,xcidx));
         else
            disp('zero len, xc=0');
         end
         
         if noiseidx==1,
            fprintf('ncount: %d %d %d %d\n',ncount(2:end));
            fprintf('(valid att conds = %d/%d [%d %d %d %d %d',...
                    length(goodatt),attcount-1,nondampedcount);
            fprintf('/%d])\n',spacelims);
         end
         % pairwise prediction test... can you find an attention
         % pair that shows more/less modulation?
         pidx=0;
         for p1=2:attcount-1,
            for p2=(p1+1):attcount,
               pidx=pidx+1;
               aattidx=find(sum(~isnan(bresp(:,end,[p1 p2])),3));
               xcp(noiseidx,1,xcidx,pidx)=xcov(rnoatt(aattidx),...
                                               ract(aattidx),0,'coeff');
               xcp(noiseidx,2,xcidx,pidx)=xcov(rdcatt(aattidx),...
                                               ract(aattidx),0,'coeff');
               xcp(noiseidx,3,xcidx,pidx)=xcov(rglobalatt(aattidx),...
                                               ract(aattidx),0,'coeff');
               xcp(noiseidx,4,xcidx,pidx)=xcov(rfullatt(aattidx),...
                                               ract(aattidx),0,'coeff');
            end
         end
         if sum(isnan(xc(noiseidx,:,xcidx))),
            keyboard
         end
         if xcidx==1,
            % save kernels
            SR(:,:,:,noiseidx,:)=lSR;
            Hpc(:,:,:,noiseidx,:)=lHpc;
         end
      end
   end
   
   predxc=cat(1,xc(1,:,:),mean(xc(2:end,:,:),1),std(xc(2:end,:,:),0,1));
   
   pxc=zeros(1,size(xc,2),size(xc,3));
   for xx=1:size(xc,2)*size(xc,3),
      nn=length(find(xc(2:end,xx)>=xc(1,xx)));
      pxc(xx)=(nn+1)./(noisecount+1);
   end
   % test for pairwise differences
   pxcp=zeros(1,size(xcp,2),size(xcp,3),size(xcp,4));
   for xx=1:size(xcp,2)*size(xcp,3)*size(xcp,4),
      nn=length(find(xcp(2:end,xx)>=xcp(1,xx)));
      pxcp(xx)=(nn+1)./(noisecount+1);
   end
   
   fprintf('XC:    %.3f %.3f %.3f %.3f\n',xc(1,:));
   fprintf('Noise: %.3f %.3f %.3f %.3f\n',...
           mean(xc(2:end,:),1));
   fprintf('PXC:   %.3f %.3f %.3f %.3f\n',pxc);
end

%targsim=zeros(attdocount,spacelim,2);
%attcc=zeros(spacelim,attuse,3);
%p=0;
pstrf=ones(spacelims,respcount,attcount);

clear x1 x2 r2 r2std resp stim eresp presp
clear data rbinset y1 y2 ystd *idx tSA timage sSAfull

