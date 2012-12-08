% cellxc2.m
%
% completes reverse correlation and stores results in these
% pre-designated matrices:
% SR = raw STA (cross correlation between signal and response)
% mH = mean decorrelated kernel (averaged over resamples)
% eH = std of decorrelated kernels across resamples
%
% (written as script to conserve memory, designed to be called by
% cellxcqueue.m. numerous parameters required) 
%
% General flow:
% 1. figure out what outputs to calucate
%    XC parameters needed: dosSA2, maxlag, boundaryfmt, stimfilter
%                          starttime,stoptime
%    normalization parms:  minegs/neigs, dotSAinvert
% 3. load appropriate stimulus/response data
%    loop to load in pieces
%    initialize output variables on first stim/loop
% 4. filter stimulus as wanted
% 5. call movxc and append results
% 6. normalize

VCELLXC=2;

% basic stim parms. should be set early
if ~exist('stimfiles'),
   disp('no stimfile selected. can''t run cellxc.m');
   keyboard
end

% parms for XC / decorr
if ~exist('maxlag'),
   %maxlag=[-10 15];
   maxlag=[0 0];
end
if ~exist('boundary'),
   boundary='zero';
end
if exist('neigs'),
   mineigs=[];
elseif exist('mineigs'),
   neigs=[];
else
   mineigs=[1 10.^(-1:-0.5:-4)];
   neigs=[];
end
% movstep: size of movie segments to send to movxc -- want this to
% be as large as possible without straining computer's memory
if ~exist('movstep'),
   movstep=50000;
end
movstep=50000;

filecount=length(respfiles);
attcount=size(times(1).start,2);  % ignore attcount>1 for the time being
attidx=1;

fidx=1;
firstseg=1;
while fidx<=filecount,
   
   % load files til nearing movstep
   curlen=0;
   stim=[];
   resp=[];
   while fidx<=filecount & ...
         (curlen==0 | stoptimes(fidx)-starttimes(fidx)+1+curlen<movstep),
      
      % load next stim file
      tstimloadparms=stimloadparms;
      if strcmp(stimloadcmd,'loadimfile') & ...
            length(stimloadparms)>0 & stimloadparms{1}==0,
         tstimloadparms{1}=round(stimloadparms{3} .* ...
                                 stimcrfs(fidx)./stimwindowcrf);
      end
      stimstart=max([starttimes(fidx,attidx)-maxlag(2) 1]);
      tstim=feval(stimloadcmd,stimfiles{fidx},stimstart,...
                     stoptimes(fidx,attidx)-maxlag(1),tstimloadparms{:});
      if ~isempty(stimfiltercmd),
         tstim=feval(stimfiltercmd,tstim,stimfilterparms{:});
      end
      
      % reshape to space X time if necessary
      iconside=size(tstim);
      bslen=iconside(end);
      iconside=iconside(1:(end-1));
      if length(iconside)>1,
         % reshape all spatial dims into one
         tstim=reshape(tstim,prod(iconside),bslen);
      end
      
      % take transpose to put time in rows and append to stim
      stim=cat(1,stim,tstim');
      
      
      % load next respfile
      disp(sprintf('Loading response: %s...',respfiles{fidx}));
      % resp=resploadatt(respfile,respvarname);
      % resp is time X space(phase) X attcode
      tresp=feval(resploadcmd,respfiles{fidx},resploadparms{:});
      
      % filter response (not used yet)
      if ~isempty(respfiltercmd),
         tresp=feval(respfiltercmd,tresp,respfilterparms{:});
      end
      
      % isolate just a single attentional set per run... reduces load
      tresp=tresp(:,:,attidx);
      
      % trim response to appropriate time range
      tresp=tresp(starttimes(fidx,attidx):stoptimes(fidx,attidx),:);
      
      rsize=size(tresp);
      tresp=tresp(:,:); % reshape to 2D matrix
      respcount=size(tresp,2);
      tresp(1:maxlag(2),:)=nan;
      lagstart=starttimes(fidx,attidx)-stimstart;
      lagend=size(tstim,2)-size(tresp,1)-lagstart;
      
      resp=cat(1,resp,ones(lagstart,respcount)*nan,tresp,...
               ones(lagend,respcount)*nan);
      
      % done, keep track
      curlen=curlen+stoptimes(fidx,attidx)-starttimes(fidx,attidx)+1;
      
      % update queue if active
      if exist('BATQUEUEID') & BATQUEUEID>0,
         dbsetqueue(BATQUEUEID,fidx);
      end
         
      fidx=fidx+1;
   end
   
   spacecount=size(stim,2);
   rsize=size(resp);
   
   meansub=1;
   if ~meansub,
      disp('subtracting mean resp!!!');
      for ii=1:size(resp,2),
         resp(:,ii)=resp(:,ii)-nanmean(resp(:,ii));
      end
   end
   
   %
   % define resampling regimes here
   %
   % option to branch according to resampfmt
   if resampcount==1 | length(resp)==0,
      % bootstrap resampling
      rstartidx=starttimes(fidx,attidx);
      rendidx=stoptimes(fidx,attidx);
   elseif resampfmt==1,
      [rstartidx,rendidx]=resampsegs(resp,resampcount);
      
   elseif resampfmt==2,
      % shuffled response resampling
      rstartidx=starttimes(fidx,attidx);
      rendidx=stoptimes(fidx,attidx);
      resp=resampnoise(resp,resampcount);
      rsize=[rsize,size(resp,3)];
      resp=resp(:,:); % reshape to 2D matrix
   end
   
   if firstseg,
      disp('initializing matrices');
      if diff(maxlag)==0 | rsize(2)==1,
         corrmtxcount=1;
         singlesSA=1;
      else
         corrmtxcount=rsize(2);
         singlesSA=0;
      end
      
      % first time, set kernel & ac matrices to zero
      SR=zeros(spacecount,diff(maxlag)+1,respcount,resampcount);
      n=zeros(respcount,resampcount);
      mS=zeros(spacecount,respcount,resampcount);
      mR=zeros(respcount,resampcount);
      tSA=zeros(diff(maxlag)*2+1,corrmtxcount,resampcount);
      sSA2=zeros(spacecount,spacecount,corrmtxcount,resampcount);
      
      % global stimulus temporal autocorr. works better than resamped?
      tSA1=zeros(diff(maxlag)*2+1,corrmtxcount,resampcount);
      nSA1=0;
      firstseg=0;
   end
   
   for resampidx=1:resampcount,
      tr=resp;
      tr([1:rstartidx(resampidx)-1 rendidx(resampidx)+1:end],:)=nan;
      %tr(rstartidx(resampidx):rendidx(resampidx),:)=nan;
      fprintf('resamp=%d: ',resampidx);
      for respidx=1:respcount,
         rgoodidx=find(~isnan(tr(:,respidx)));
         tmR=sum(resp(rgoodidx,respidx))';
         tmS=sum(stim(rgoodidx,:))';
         tn=length(rgoodidx);
         
         tSR=zeros(spacecount,diff(maxlag)+1);
         for tt=maxlag(1):maxlag(2),
            fprintf('.');
            trg=rgoodidx;
            trg=trg(find(trg-tt>0 & trg-tt<size(stim,1)));
            tSR(:,tt-maxlag(1)+1)=(resp(trg,respidx)'*stim(trg-tt,:))';
         end
         
         if respidx<=corrmtxcount,
            if dosSA==2,
               tsSA2=stim(rgoodidx,:)'*stim(rgoodidx,:);
               fprintf('*');
            end
            for xx=1:spacecount,
               tstim=stim(rgoodidx,xx);
               ttSA=xcorr(tstim,diff(maxlag),'biased') ./ ...
                    spacecount.*length(rgoodidx);
            end
            fprintf('*\n');
         end
      end
      
      for rr=1:resampcount,
         % add means to all resamp channels. this is to
         % make it so that everything is centered around
         % the same DC!
         if resampcount==1 | rr~=resampidx,
            % otherwise add outputs to running total
            SR(:,:,respidx,rr)=SR(:,:,respidx,rr)+tSR;
            mS(:,respidx,rr)=mS(:,respidx,rr)+tmS;
            mR(respidx,rr)=mR(respidx,rr)+tmR;
            n(respidx,rr)=n(respidx,rr)+tn;
            tSA(:,:,rr)=tSA(:,:,rr)+ttSA;
            
            if dosSA==2 & respidx<=corrmtxcount,
               sSA2(:,:,respidx,rr)=sSA2(:,:,respidx,rr)+tsSA2;
            end
         end
      end
      
      % update queue if active
      if exist('BATQUEUEID') & BATQUEUEID>0,
         dbsetqueue(BATQUEUEID,resampidx*10);
      end
   end
   
   for xx=1:spacecount,
      tstim=stim(:,xx)-mean(stim(:,xx));
      tSA1=tSA1+repmat(xcorr(tstim,diff(maxlag),'biased')./ ...
                       spacecount.*size(stim,1),...
                       [1 corrmtxcount resampcount]);
   end
   nSA1=nSA1+size(stim,1);
end

% zeroth order normalization. ie, divide sums by the number of
% samples to get appropriate means
[SR,mS,mR,tSA,sSA2]=normalize0(SR,n,mS,mR,tSA,sSA2,meansub);
tSA1=tSA1./nSA1;

if dosSA==3,
   load ([stimfiles{fidx},'.AC.mat'],'sSA2');
elseif 0,
   disp('taking SINGLE sSA2!');
   % take average sSA2 to speed things up and possibly reduce
   % noise in kernel estimates... hell, it worked in V4
   sSA2=mean(sSA2,4);
end
%[H,neigs]=normalize(SR,sSA2,tSA,mineigs,neigs);
%[H,lambda]=normalizereg(SR,sSA2,tSA,sfscount);

% 8/7/02: use single tSA for all resamples
[H,lambda]=normalizereg(SR,sSA2,tSA1,sfscount,[],1); % topSR=1
%sfscount=10;
%[H,lambda]=normalizeshr(SR,sSA2,tSA1,sfscount);
%[H,neigs]=normalize(SR,sSA2,tSA1,mineigs,neigs);

irS=squeeze(std(SR,1,1));
irH=squeeze(std(H,1,1));

allt=std(reshape(permute(H,[1 2 4 5 3]),spacecount*(diff(maxlag)+1)* ...
                 respcount*resampcount,sfscount));
alls=std(reshape(permute(H,[1 3 4 5 2]),spacecount*sfscount* ...
                 respcount*resampcount,(diff(maxlag)+1)));
tH=squeeze(std(reshape(permute(H,[1 4 5 2 3]),spacecount* ...
                       respcount*resampcount,(diff(maxlag)+1),sfscount)));


mH=mean(H,5);
eH=std(H,1,5) .* sqrt((resampcount-1)/resampcount);

mSR=mean(SR,4);
eSR=std(SR,1,4).* sqrt((resampcount-1)/resampcount);

% take mean stim/response across resamples
mSall=mean(mS,3);
mRall=mean(mR,2);
ntot=n;

if DOCELLFIT2,
   cellfit2;
end

disp('cellxc2 done?');
return


% 4/1/02
% new ALGORITHM to optimize processing speed and save memory:
% figure out how many attentional states there are.
% for each att state:
%   for each file
%      load response data for current att state
%      figure out resample segments over length of all movies
%      for each resample segment (if resamples longer than movstep):
%         for each movseg
%            accumulate SR/SA info
%      assemble resampled kernels (or temporally indep ones???)
%      if maxlag=[0 0], reshape phase/space dim to "latency" dim
%      decorrelate
%      save mean/variance kernel for each decorr specification and
%         clear the individual resampled kernels
% then do preds -- cellfit.m, cellpred.m

% each dimension of the response is special:
%    time X phase(space) X attention(temporal subset) [X resample]

% set up matrices to store decorrelated kernels
mH=[];    % mean across resamples
smH=[];   % mean across smoothed resamples
eH=[];    % std across smoothed resamples
mSR=[];   % mean STA across resamples
eSR=[];   % std of STA across resamples
ntot=[];
mSall=[]; % mean stimulus - for subtracting in preds.
mRall=[]; % mean response



% for each att state:
for attidx=1:attcount,
   fcount=max(find(starttimes(:,attidx)<stoptimes(:,attidx)));
   for fidx=1:fcount,
      
      
      %
      % define resampling regimes here
      %
      % option to branch according to resampfmt
      if resampcount==1 | length(resp)==0,
         % bootstrap resampling
         rstartidx=starttimes(fidx,attidx);
         rendidx=stoptimes(fidx,attidx);
      elseif resampfmt==1,
         [rstartidx,rendidx]=resampsegs(resp,resampcount);
         if strcmp(boundary,'zero'),
            rstartidx=rstartidx-maxlag(2);
            rendidx=rendidx-maxlag(1);
            rstartidx(find(rstartidx<1))=1;
            if 0 & USEBSTIM,
               % disabled svd 8/8/02... doesn't help?
               rendidx(find(rendidx>size(bstim,1)))=stoptimes(fidx,attidx);
            else
               rendidx(find(rendidx>stoptimes(fidx,attidx)))=...
                   stoptimes(fidx,attidx);
            end
         end
         
      elseif resampfmt==2,
         % shuffled response resampling
         rstartidx=starttimes(fidx,attidx);
         rendidx=stoptimes(fidx,attidx);
         resp=resampnoise(resp,resampcount);
         rsize=[rsize,size(resp,3)];
         resp=resp(:,:); % reshape to 2D matrix
      end
      
      % extend resp with nan's to allow for extra normalizat?????
      resp(size(resp,1)+1:rendidx(end),:)=nan;
      resptot=size(resp,2);
      
      % for PFTH RC, only do a single sSA for all response channels,
      % otherwise do a separate one for each response channel
      if diff(maxlag)==0 & rsize(2)>1,
         corrmtxcount=1;
         singlesSA=1;
      else
         corrmtxcount=resptot;
         singlesSA=0;
      end
      
      for resampidx=1:length(rstartidx),
         % figure out parms for stimulus segmenting
         curstimidx=rstartidx(resampidx)-1;
         len=rendidx(resampidx)-rstartidx(resampidx)+1; 
         movsegcount=ceil(len/movstep);
         
         % loop through stimulus segments and accumulate their
         % contributions to the kernel and autocorrelation estimates.
         for movsegidx=1:movsegcount,
            startidx=curstimidx+1;
            stopidx=min([curstimidx+movstep rendidx(resampidx)]);
            
            if exist('USEBSTIM') & USEBSTIM,
               stim=bstim(startidx:stopidx,:);
            else
               % if scale_pix not defined, use crf parms
               tstimloadparms=stimloadparms;
               if strcmp(stimloadcmd,'loadimfile') & ...
                     length(stimloadparms)>0 & stimloadparms{1}==0,
                  tstimloadparms{1}=stimloadparms{3} .* ...
                      stimcrfs(fidx)./stimwindowcrf;
               end
               stim=feval(stimloadcmd,stimfiles{fidx},startidx,stopidx,...
                          tstimloadparms{:});
               
               % filter stimulus segment if selected
               if ~isempty(stimfiltercmd),
                  stim=feval(stimfiltercmd,stim,stimfilterparms{:});
               end
               
               % reshape to space X time if necessary
               iconside=size(stim);
               iconside=iconside(1:(end-1));
               if length(iconside)>1,
                  % reshape all spatial dims into one
                  stim=reshape(stim,prod(iconside),stopidx-startidx+1);
               end
               stim=stim'; % take transpose to put time in rows
               spacecount=size(stim,2);
            end
            
            if stopidx-startidx==0,
               tSR=zeros(spacecount,diff(maxlag)+1,resptot);
               tn=zeros(resptot,1);
               ttSA=zeros(diff(maxlag)*2+1,corrmtxcount);
               tsSA2=zeros(spacecount,spacecount,corrmtxcount);
               
            elseif dosSA==2,
               % meansub is forced to be 0 here. real meansub is
               % completed in normalize0
               [tSR,tn,tmS,tmR,ttSA,tsSA2]=...
                   movxc(stim,resp(startidx:stopidx,:),maxlag,boundary,...
                         singlesSA,0);
            else
               % if dosSA==3, don't calc sSA2, just load it later
               [tSR,tn,tmS,tmR,ttSA]=movxc(stim,resp(startidx:stopidx,:),...
                                           maxlag,boundary,singlesSA,0);
               tsSA2=[];
            end
            
            if resampfmt==1,
               if movsegidx==1 & resampidx==1 & fidx==1,
                  % first time, set kernel & ac matrices to zero
                  clear SR n tSA sSA2
                  SR=zeros(spacecount,size(tSR,2),respcount,resampcount);
                  n=zeros(respcount,resampcount);
                  mS=zeros(spacecount,respcount,resampcount);
                  mR=zeros(respcount,resampcount);
                  tSA=zeros(size(ttSA,1),corrmtxcount,resampcount);
                  %sSA1=zeros(spacecount,spacecount,corrmtxcount);
                  %nSA1=0;
                  sSA2=zeros(spacecount,spacecount,corrmtxcount,resampcount);
               end
               for rr=1:resampcount,
                  % add means to all resamp channels. this is to
                  % make it so that everything is centered around
                  % the same DC!
                  if resampcount==1 | rr~=resampidx,
                     % otherwise add outputs to running total
                     SR(:,:,:,rr)=SR(:,:,:,rr)+tSR;
                     mS(:,:,rr)=mS(:,:,rr)+tmS;
                     mR(:,rr)=mR(:,rr)+tmR;
                     n(:,rr)=n(:,rr)+tn;
                     tSA(:,:,rr)=tSA(:,:,rr)+ttSA;
                     %sSA1=sSA1+tsSA2;
                     %nSA1=nSA1+tn;
                     
                     if dosSA==2,
                        sSA2(:,:,:,rr)=sSA2(:,:,:,rr)+tsSA2;
                     end
                  end
               end
            elseif resampfmt==2,
               if movsegidx==1 & resampidx==1 & fidx==1,
                  % first time, set kernel & ac matrices to zero
                  clear SR n tSA sSA2
                  SR=zeros(spacecount,size(tSR,2),respcount,resampcount);
                  n=zeros(respcount,resampcount);
                  mS=zeros(spacecount,respcount,resampcount);
                  mR=zeros(respcount,resampcount);
                  tSA=zeros(size(ttSA,1),corrmtxcount);
                  sSA2=zeros(spacecount,spacecount,corrmtxcount);
               end
               % otherwise add outputs to running total
               SR=SR+reshape(tSR,spacecount,size(tSR,2),respcount,resampcount);
               n=n+reshape(tn,respcount,resampcount);
               mS=mS+reshape(tmS,spacecount,respcount,resampcount);
               mR=mR+reshape(tmR,respcount,resampcount);
               tSA=tSA+reshape(ttSA,size(ttSA,1),corrmtxcount);
               if dosSA==2,
                  sSA2=sSA2+reshape(tsSA2,spacecount,spacecount,corrmtxcount);
               end
            end
            
            curstimidx=curstimidx+movstep;
         end
         
         fprintf('Completed att state %d, file set %d, resamp seg %d.\n',...
                 attidx,fidx,resampidx);
         
         % update queue if active
         if exist('BATQUEUEID') & BATQUEUEID>0,
            dbsetqueue(BATQUEUEID,(attidx-1)*1000+(fidx-1)*100+resampidx);
         end
      end
      
      % tSA compute hack: tSA1 is same for all resample
      % conditions. this seems to give better predictions.
      if movsegidx==1 & fidx==1,
         tSA1=zeros(size(ttSA,1),corrmtxcount,resampcount);
         nSA1=0;
      end
      %userange=starttimes(fidx,attidx):stoptimes(fidx,attidx);
      userange=rstartidx(1):rendidx(end);
      for xx=1:spacecount,
         if meansub,
            tstim=bstim(userange,xx)-mean(bstim(userange,xx));
         else
            tstim=bstim(userange,xx);
         end
         tSA1=tSA1+repmat(xcorr(tstim,diff(maxlag),'biased')./ ...
                        spacecount.*length(userange),...
                        [1 corrmtxcount resampcount]);
      end
      nSA1=nSA1+length(userange);
   end
   
   % divide each thing by the number of samples. 0th-order
   % normalization. also subtract the mean ... here rather than
   % each separate movxc piece
   % function [SR,mS,mR,tSA,sSA2]=normalize0(SR,n,mS,mR,tSA,sSA2)
   [SR,mS,mR,tSA,sSA2]=normalize0(SR,n,mS,mR,tSA,sSA2,meansub);
   %sSA1=sSA1./nSA1;
   tSA1=tSA1./nSA1;
   
   % then do 2nd order normalization
   % (for PFTH analysis, reshape kernels to be more convenient)
   if ~dosSA,
      [H,lambda]=normalizetime(SR,tSA,sfscount);
      neigs=sfscount;
   elseif dosSA==1,
      [H,neigs,lambda]=normalize1(SR,sSA1,tSA,mineigs,neigs);
      
   elseif singlesSA,
      % this is the case where resampled kernels are noise. want to
      % do the decorrelation in the eigenvector domain of the full
      % stimulus but corrected for correlation in the specific
      % attentionally restricted stimulus sequence
      
      SEPDOMS=1;
      if ~SEPDOMS,
         if attidx==1,
            % average AC matrix of full attentional state over resamples
            sSA0=mean(sSA2(:,:,1,:),4);
            [u0,s0,v0]=svd(sSA0);

            ds0=1./diag(s0);
            ds0(find(ds0>(ds0(1)*10^10)))=0;
         end
         
         SR=reshape(SR,spacecount,respcount,1,resampcount);
         tSR=SR;
         
         % decorrelate SR wrt sSA0 correlations. should be modest
         % and relatively unchanged
         SR=reshape(SR,spacecount,respcount,1,resampcount);
         tSR=zeros(size(SR));
         for resampidx=1:corrmtxcount*resampcount,
            tSR(:,:,1,resampidx)=(u0*diag(ds0)*u0'*sSA2(:,:,1,resampidx)) ...
                * tSR(:,:,1,resampidx);
         end
      else
         SR=reshape(SR,spacecount,respcount,1,resampcount);
         tSR=SR;
         sSA0=sSA2;
      end
      
      % NOTICE: normalize tSR with sSA0--not SR with sSA2!
      %[H,neigs]=normalize(tSR,sSA0,[],mineigs,neigs);
      %[H,lambda]=normalizereg(tSR,sSA0,[],sfscount,[],0);
      %[H,lambda]=normalizereg(tSR,sSA0,[],sfscount,[],1); % topSR=1
      [H,lambda]=normalizereg(tSR,sSA0,[],sfscount,[],0); % topSR=0
   else
      if dosSA==3,
         load ([stimfiles{fidx},'.AC.mat'],'sSA2');
      elseif 0,
         disp('taking SINGLE sSA2!');
         % take average sSA2 to speed things up and possibly reduce
         % noise in kernel estimates... hell, it worked in V4
         sSA2=mean(sSA2,4);
      end
      %[H,neigs]=normalize(SR,sSA2,tSA,mineigs,neigs);
      %[H,lambda]=normalizereg(SR,sSA2,tSA,sfscount);
      
      % 8/7/02: use single tSA for all resamples
      [H,lambda]=normalizereg(SR,sSA2,tSA1,sfscount,[],1); % topSR=1
      %sfscount=10;
      %[H,lambda]=normalizeshr(SR,sSA2,tSA1,sfscount);
      %[H,neigs]=normalize(SR,sSA2,tSA1,mineigs,neigs);
   end
   
   irS=squeeze(std(SR,1,1));
   irH=squeeze(std(H,1,1));
   
   DOSMOOTHORSF=0;
   %sH=H;
   %keyboard
   if DOSMOOTHORSF & (strcmp(kernfmt,'fft') | strcmp(kernfmt,'pfft')),
      fprintf('OR/SF smoothing H');
      
      for ridx=1:resampcount,
         fprintf(' resamp=%d...',ridx);
         if strcmp(kernfmt,'fft') ,
            phasecount=4;
         else
            phasecount=1;
         end
         sfbincount=sqrt(spacecount/phasecount*2)/2;
         orbincount=round(12*sfbincount/8);
         
         H(:,-maxlag(1)+1:end,:,:,ridx)=...
             kernsmoothorsf(H(:,-maxlag(1)+1:end,:,:,ridx),orbincount,...
                            sfbincount,phasecount,7,0.3,0.3);  
         % .7 both was best for batch 23

         %sH(:,-maxlag(1)+1:end,:,:,ridx)=H(:,-maxlag(1)+1:end,:,:,ridx);
         %sH(:,-maxlag(1)+1:end,:,:,ridx)=...
         %    kernsmoothorsf(sH(:,-maxlag(1)+1:end,:,:,ridx),8,4,phasecount,11);
         if exist('BATQUEUEID') & BATQUEUEID>0,
            dbsetqueue(BATQUEUEID,attcount*1000+ridx*100);
         end
      end
   end
   keyboard
   
   % calculate means and variance across resampled kernels
   % append to collection of kernels for different attentional states
   if resampcount>1 & resampfmt==1,  
      % might want to add other logic here for
      % handling things other than bootstrapping
      % for time being ... just compute standard error with
      % bootstrapping formula
      mH=cat(5,mH,mean(H,5));
      eH=cat(5, eH, std(H,1,5).*sqrt((resampcount-1)/resampcount));
      %smH=cat(5,smH,mean(sH,5));
      %eH=cat(5,eH,std(sH,1,5)).* sqrt((resampcount-1)/resampcount);
      mSR=cat(4,mSR,mean(SR,4));
      eSR=cat(4,eSR,std(SR,1,4)).* sqrt((resampcount-1)/resampcount);
   elseif resampcount>1 & resampfmt==2,  
      % might want to add other logic here for
      % handling things other than bootstrapping
      % for time being ... just compute standard error with
      % bootstrapping formula
      mH=cat(5,mH,H(:,:,:,:,1));
      eH=cat(5,eH,std(H(:,:,:,:,2:end),1,5));
      %smH=cat(5,smH,sH(:,:,:,:,1));
      %eH=cat(5,eH,std(sH(:,:,:,:,2:end),1,5));
      mSR=cat(4,mSR,mean(SR,4));
      eSR=cat(4,eSR,std(SR,1,4)).* sqrt((resampcount-1)/resampcount);
      tmH=mean(H(:,:,:,:,2:end),5);
   else
      mH=cat(5,mH,H);
      %smH=cat(5,smH,sH);
      eH=cat(5,eH,0);
      mSR=cat(4,mSR,SR);
      eSR=cat(4,eSR,0);
   end
   
   % take mean stim/response across resamples
   mSall=cat(3,mSall,mean(mS,3));
   mRall=cat(2,mRall,mean(mR,2));
   ntot=cat(3,ntot,n);
   
   if DOCELLFIT2,
      cellfit2;
   end
   
   clear H sH
   
   %figure(1);
   %showkern(mH(:,(1-maxlag(1)):end,round(linspace(1,length(neigs),6))),kernfmt,iconside);
   %showkern(mH(:,(1-maxlag(1)):end,:),kernfmt,iconside);
   %showkern(squeeze(mH(:,(1-maxlag(1)):end,1,1,:)./sH(:,(1-maxlag(1)):end,1,1,:)),kernfmt);
   %showkern(SR(:,(1-maxlag(1)):end,1:5),kernfmt,iconside);
   %drawnow;
end

if spacecount>10,
   % clear unwanted space wasters
   clear sH H ttSA tsSA2 tSR
   clear tSA sSA1 sSA2 sSA0 SR
end

clear stim resp bstim

