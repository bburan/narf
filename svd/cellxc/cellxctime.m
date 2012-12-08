% cellxc.m
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
   movstep=10000;
end
if ~exist('DOCELLFIT2'),
   DOCELLFIT2=0;
end

filecount=length(respfiles);
attcount=size(times(1).start,2);

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
%
% 11/21/02: pruned out attention and PFTH functionality.. want
% streamlined, optimized (bug free!) estimation procedure
%

% each dimension of the response is special:
%    time X phase(space) X attention(temporal subset) [X resample]

if ~exist('attidx','var'),
   attidx=1;
end
fcount=max(find(starttimes(:,attidx)<stoptimes(:,attidx)));
for fidx=1:fcount,
      
   % check to see if stim and response exist
   if ~checkbic & not(exist(respfiles{fidx},'file')),
      disp([respfiles{fidx},' missing!!']);
      return
   end
   if ~checkbic & not(exist(stimfiles{fidx},'file')),
      disp([stimfiles{fidx},' missing!!']);
      return
   end
   
   %
   % load the response data
   %
   disp(sprintf('Loading response: %s...',respfiles{fidx}));

   % resp=resploadatt(respfile,respvarname);
   % resp is time X space(phase) X attcode
   resp=feval(resploadcmd,respfiles{fidx},resploadparms{:});
      
   % filter response (not used yet)
   if ~isempty(respfiltercmd),
      resp=feval(respfiltercmd,resp,respfilterparms{:});
   end
      
   % isolate just a single attentional set per run... reduces load
   bresp=resp;
   resp=resp(:,:,attidx);
   
   % trim response to appropriate time range
   resp=resp(1:stoptimes(fidx,attidx),:);
   
   rsize=size(resp);
   resp=resp(:,:); % reshape to 2D matrix
   respcount=size(resp,2);
      
   % LOAD WHOLE STIMULUS -- speeds things up and old-style
   % mean sub works for short movies. very hacky. need to
   % document better. how?
   USEBSTIM=1;
   if USEBSTIM,
      tstimloadparms=stimloadparms;
      if strcmp(stimloadcmd,'loadimfile') & ...
            length(stimloadparms)>0 & stimloadparms{1}==0,
         tstimloadparms{1}=round(stimloadparms{3} .* ...
                                 stimcrfs(fidx)./stimwindowcrf);
      end
      bstim=feval(stimloadcmd,stimfiles{fidx},1,...
                  stoptimes(fidx,attidx)-maxlag(1),tstimloadparms{:});
      if ~isempty(stimfiltercmd),
         bstim=feval(stimfiltercmd,bstim,stimfilterparms{:});
      end
      
      % reshape to space X time if necessary
      iconside=size(bstim);
      bslen=iconside(end);
      iconside=iconside(1:(end-1));
      if length(iconside)>1,
         % reshape all spatial dims into one
         bstim=reshape(bstim,prod(iconside),bslen);
      end
      bstim=bstim'; % take transpose to put time in rows
      spacecount=size(bstim,2);
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
   
   % extend resp with nan's to allow for extra normalization?????
   resp(size(resp,1)+1:rendidx(end),:)=nan;
   resptot=size(resp,2);
   
   corrmtxcount=resptot;
   singlesSA=0;
      
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
         
         curstimidx=curstimidx+movstep;
      end
      
      fprintf('Completed file set %d, resamp seg %d.\n',fidx,resampidx);
      
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
   
else
   if dosSA==3,
      load ([stimfiles{fidx},'.AC.mat'],'sSA2');
   end
   
   %[H,neigs]=normalize(SR,sSA2,tSA,mineigs,neigs);
   %[H,lambda]=normalizereg(SR,sSA2,tSA,sfscount);
   
   % 8/7/02: use single tSA for all resamples
   %keyboard
   [H,lambda]=normalizereg(SR,sSA2,tSA1,sfscount);
   %sfscount=10;
   %[H,lambda]=normalizeshr(SR,sSA2,tSA1,sfscount);
   %[H,neigs]=normalize(SR,sSA2,tSA1,mineigs,neigs);
end

irS=squeeze(std(SR,1,1));
irH=squeeze(std(H,1,1));

DOSMOOTHORSF=1;
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
                         sfbincount,phasecount,7,0.6,0.6);
      %sH(:,-maxlag(1)+1:end,:,:,ridx)=H(:,-maxlag(1)+1:end,:,:,ridx);
      %sH(:,-maxlag(1)+1:end,:,:,ridx)=...
      %    kernsmoothorsf(sH(:,-maxlag(1)+1:end,:,:,ridx),8,4,phasecount,11);
      if exist('BATQUEUEID') & BATQUEUEID>0,
         dbsetqueue(BATQUEUEID,attcount*1000+ridx*100);
      end
   end
end

% calculate means and variance across resampled kernels
% append to collection of kernels for different attentional states
if resampcount>1,  
   % might want to add other logic here for
   % handling things other than bootstrapping
   % for time being ... just compute standard error with
   % bootstrapping formula
   mH=mean(H,5);
   eH=std(H,1,5).*sqrt((resampcount-1)/resampcount);
   mSR=mean(SR,4);
   eSR=std(SR,1,4).* sqrt((resampcount-1)/resampcount);
else
   mH=H;
   eH=0;
   mSR=cat(4,mSR,SR);
   eSR=cat(4,eSR,0);
end

% take mean stim/response across resamples
mSall=mean(mS,3);
mRall=mean(mR,2);
ntot=n;

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

if spacecount>10,
   % clear unwanted space wasters
   clear sH H ttSA tsSA2 tSR
   clear tSA sSA1 sSA2 sSA0 SR
end

clear stim resp bstim

