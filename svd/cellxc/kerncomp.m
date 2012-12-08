% kerncomp.m : compare kernels from different attentional states to
% see how similar they are

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

% attuse: hacky parameter telling how many attentional conditions
% there are (4 targets plus "all/-1")
attuse=5;

% noisecount already defined... number of fake "attentional" states
attcount=attuse+noisecount;
if strcmp(respfiltercmd,''),
   resploadparms{5}=noisecount;
else
   respfilterparms{3}=noisecount;
end

% for each att state:
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
   rsize=size(bresp);
   
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
      starttimes(fidx,attidx)=rvalid{fidx,attidx}(1);
      stoptimes(fidx,attidx)=rvalid{fidx,attidx}(resplens(fidx,attidx));
   end
   
   brmean=nanmean(bresp(:,:,1));
   
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
   
   
   for attidx=1:attcount,
      resp=bresp(:,:,attidx);
      
      % trim response to appropriate time range
      resp=resp(1:stoptimes(fidx,attidx),:);
      
      rsize=size(resp);
      resp=resp(:,:); % reshape to 2D matrix
      respcount=size(resp,2);
      
      if respfmtcode==0,
         tbincount=diff(maxlag)+1;
      else
         tbincount=respcount;
      end
      %
      % define resampling regimes here
      %
      % option to branch according to resampfmt
      if resampcount>1,
         [rstartidx,rendidx]=resampsegs(resp,resampcount);
         if strcmp(boundary,'zero'),
            rstartidx=rstartidx-maxlag(2);
            rendidx=rendidx-maxlag(1);
            rstartidx(find(rstartidx<1))=1;
            rendidx(find(rendidx>stoptimes(fidx,attidx)))=...
                stoptimes(fidx,attidx);
         end
      else
         rstartidx=1;
         rendidx=size(resp,1);
      end
      
      % for PFTH RC, only do a single sSA for all response channels,
      % otherwise do a separate one for each response channel
      if diff(maxlag)==0 & rsize(2)>1,
         corrmtxcount=1;
         singlesSA=1;
      else
         corrmtxcount=respcount;
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
            
            if PRELOADMOV==1,
               fprintf('Using stim %d-%d...\n',startidx,stopidx);
               stim=bstim(startidx:stopidx,:);
            else
               stim=feval(stimloadcmd,tstimfile,startidx,stopidx,...
                          stimloadparms{:});
               
               % filter stimulus segment if selected
               if ~isempty(stimfiltercmd),
                  stim=feval(stimfiltercmd,stim,stimfilterparms{:});
               end
               % reshape to space X time if necessary
               iconside=size(stim);
               iconside=iconside(1:(end-1));
               if length(iconside)>=2,
                  % reshape all spatial dims into one
                  stim=reshape(stim,prod(iconside),stopidx-startidx+1);
               end
               stim=stim'; % take transpose to put time in rows
               spacecount=size(stim,2);
            end
            
            if stopidx-startidx==0,
               tSR=zeros(spacecount,diff(maxlag)+1,respcount);
               tn=zeros(respcount,1);
               tsSA1=zeros(spacecount,corrmtxcount);
               ttSA=zeros(diff(maxlag)*2+1,corrmtxcount);
               tsSA2=zeros(spacecount,spacecount,corrmtxcount);
               
            elseif dosSA>1,
               [tSR,tn,tmS,tmR,ttSA,tsSA2]=...
                   movxc(stim,resp(startidx:stopidx,:),maxlag,boundary,...
                         singlesSA,0);
            else
               [tSR,tn,tmS,tmR,ttSA]=...
                   movxc(stim,resp(startidx:stopidx,:),maxlag,boundary,...
                         singlesSA,0);
               tsSA2=[];
            end
            
            if movsegidx==1 & resampidx==1,
               % first time, set kernel & ac matrices to zero
               clear SR n sSA1 tSA sSA2
               SR=zeros(spacecount,size(tSR,2),respcount,resampcount);
               n=zeros(respcount,resampcount);
               mS=zeros(spacecount,respcount,resampcount);
               mR=zeros(respcount,resampcount);
               tSA=zeros(size(ttSA,1),corrmtxcount,resampcount);
               sSA2=zeros(spacecount,spacecount,corrmtxcount,resampcount);
            end
            
            for rr=1:resampcount,
               if resampcount==1 | rr~=resampidx,
                  % otherwise add outputs to running total
                  SR(:,:,:,rr)=SR(:,:,:,rr)+tSR;
                  n(:,rr)=n(:,rr)+tn;
                  mS(:,:,rr)=mS(:,:,rr)+tmS;
                  mR(:,rr)=mR(:,rr)+tmR;
                  tSA(:,:,rr)=tSA(:,:,rr)+ttSA;
                  if dosSA>1,
                     sSA2(:,:,:,rr)=sSA2(:,:,:,rr)+tsSA2;
                  end
               end
            end
            curstimidx=curstimidx+movstep;
         end
         
         fprintf('Done att%d f%d resamp%d.\n',attidx,fidx,resampidx);
         
         % update queue if active
         if exist('BATQUEUEID') & BATQUEUEID>0,
            dbsetqueue(BATQUEUEID,(attidx-1)*1000+(fidx-1)*100+resampidx);
         end
      end
      
      % RC components have been calculated for current attentional
      % state, now normalize them and convert to eigenvector
      % domain.
      if attidx==1,
         if respfmtcode==0,
            spacelim=15;  % spatial bins to save in eigS/eigH
            ktestcount=spacelim;        % number of regularization parms
         elseif DOREGNORM ,
            spacelim=spacecount-1;  % spatial bins to save in eigS/eigH
            ktestcount=30;        % number of regularization parms
         else
            spacelim=30;          % spatial bins to save in eigS/eigH
            if spacelim>spacecount,
               spacelim=spacecount-1;
            end
            ktestcount=spacelim;  % do T2 on 1 through spacelim eigs
         end
         
         DOPAIRS=1;
         if DOPAIRS,
            attpairs=[1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
            attdocount=6;
         else
            attdocount=4;
         end
         
         noisesets=floor((noisecount+attuse-1)./(attuse-1));
         if respcount>1,
            resptouse=respcount;  % don't use all respcount eigH's
            respbase=0;
         else
            resptouse=1;
            respbase=1;
         end
         T2=zeros(attdocount,noisesets,ktestcount,resptouse);
         p=ones(attdocount,ktestcount,resptouse);
         
         % these get saved for export
         eigS=zeros(spacelim,tbincount,attcount);
         eigSout=zeros(spacelim,tbincount,attcount);
         eigH=zeros(spacelim,tbincount,attcount);
         eigHout=zeros(spacelim,tbincount,attcount);
         
         % for first attentional state (all attention) create
         % eigenvector domain matrices with appropriate sizes
         eigSatt=zeros(spacelim,tbincount,resampcount,attuse-1);
         eigSoutatt=zeros(spacelim,tbincount,resampcount,attuse-1);
         eigHatt=zeros(spacelim,tbincount,resampcount,attuse-1);
         eigHoutatt=zeros(spacelim,tbincount,resampcount,attuse-1);
         
         % compute mean stimulus and response for each resampled
         % kernel then subtract mean, normalize by sample count the
         % raw STA (SR) and the autocorrelation matrix (sSA2)
         [SR,mS,mR,tSA,sSA2]=normalize0(SR,n,mS,mR,tSA,sSA2,1);
         
         % compute eigenvectors from AC matrix of entire stimulus sequence
         SA=sSA2(:,:,1,1);
         [U0,S0,v]=svd(SA);
         S0=diag(S0);
         
         uu=U0(:,1:spacelim);
         sd=S0(1:spacelim);
         
         % project STA into eigenvector domain
         for resampidx=1:resampcount,
            eigSatt(:,:,resampidx,1)=...
                   uu' * squeeze(SR(:,:,:,resampidx));
            for eigidx=1:spacelim,
               eigHatt(eigidx,:,resampidx,1)=...
                   eigSatt(eigidx,:,resampidx,1)./sd(eigidx);
            end
         end
         eigS(:,:,1)=mean(eigSatt(:,:,:,1),3);
         eigH(:,:,1)=mean(eigHatt(:,:,:,1),3);
         eigHerr=std(eigHatt(:,:,:,1),1,3).*sqrt((resampcount-1)/resampcount);
   
         snrthresh=10;
         
         snr=abs(eigH(:,:,1)./eigHerr);
         sigidx=snr>snrthresh;
         
         inidx=1;
      else
         % for attidx>1
         % accumulate appropriate in/out attentional data in
         % current set of 4 SR/AC matrices 
         % "in" kernel is current attidx, also gets added to "out"
         % kernel of other 3 complementary attentional states
         sampidx=ceil((attidx-1)/(attuse-1));
         sbase=(sampidx-1)*(attuse-1)+1;
         inidx=attidx-sbase;
         
         if inidx==1,
            SRatt=zeros([size(SR) attuse-1 2]);
            SAatt=zeros([size(sSA2) attuse-1 2]);
            natt=zeros([size(n) attuse-1 2]);
            mSatt=zeros([size(mS) attuse-1 2]);
            mRatt=zeros([size(mR) attuse-1 2]);
         end
         SRatt(:,:,:,:,inidx,1)=SR;
         SAatt(:,:,:,:,inidx,1)=sSA2;
         natt(:,:,inidx,1)=n;
         mSatt(:,:,:,inidx,1)=mS;
         mRatt(:,:,inidx,1)=mR;
         
         for outidx=[1:(inidx-1) (inidx+1):(attuse-1)],
            SRatt(:,:,:,:,outidx,2)=SRatt(:,:,:,:,outidx,2)+SR;
            SAatt(:,:,:,:,outidx,2)=SAatt(:,:,:,:,outidx,2)+sSA2;
            natt(:,:,outidx,2)=natt(:,:,outidx,2)+n;
            mSatt(:,:,:,outidx,2)=mSatt(:,:,:,outidx,2)+mS;
            mRatt(:,:,outidx,2)=mRatt(:,:,outidx,2)+mR;
         end
      end
      
      if inidx==attuse-1, % process into eigenvector domain
         
         fprintf('Converting attidx=%d-%d to PC-domain...\n',...
                 sbase+1,sbase+attuse-1);
         
         % normalize everything as needed
         for inidx=1:(attuse-1)*2,
            for resampidx=1:resampcount,
               mSatt(:,1,resampidx,inidx)=...
                   mSatt(:,1,resampidx,inidx)./natt(1,resampidx,inidx);
               for latidx=1:respcount,
                  mRatt(latidx,resampidx,inidx)=...
                      mRatt(latidx,resampidx,inidx)./...
                      natt(latidx,resampidx,inidx);
                  for tidx=1:diff(maxlag)+1,
                     SRatt(:,tidx,latidx,resampidx,inidx)=...
                         SRatt(:,tidx,latidx,resampidx,inidx)./...
                         natt(latidx,resampidx,inidx) -...
                         mSatt(:,1,resampidx,inidx).* ...
                         mRatt(latidx,resampidx,inidx);
                  end
               end
               SAatt(:,:,1,resampidx,inidx)=...
                   SAatt(:,:,1,resampidx,inidx)./natt(1,resampidx,inidx)-...
                   mSatt(:,1,resampidx,inidx)*mSatt(:,1,resampidx,inidx)';
            end
         end
         
         if respfmtcode==0 & exist('ttSAinv'),
            totkcount=prod(size(natt));
            for respidx=1:totkcount,
               SRatt(:,:,respidx)=SRatt(:,:,respidx)*ttSAinv';
               SRatt(:,:,respidx)=conv2(SRatt(:,:,respidx),...
                                        [0.2 0.6 0.2],'same');
            end
         end
         
         % project into eignevector domain of "all" atten condition
         for inidx=1:(attuse-1),
            for resampidx=1:resampcount,
               eigSatt(:,:,resampidx,inidx)=...
                   uu' * squeeze(SRatt(:,:,:,resampidx,inidx,1));
               eigSoutatt(:,:,resampidx,inidx)=...
                   uu' * squeeze(SRatt(:,:,:,resampidx,inidx,2));
               
               if 0,
                  % SVD 6/21/02: turn off 2nd-order decorr. decorr, though
                  % theoretically better, seems to make things worse?
                  B=diag(sd(:,1));
                  eigHatt(:,:,resampidx,inidx)=...
                      B^-1 * eigSatt(:,:,resampidx,inidx);
                  eigHoutatt(:,:,resampidx,inidx)=...
                      B^-1 * eigSoutatt(:,:,resampidx,inidx);
               else
                  % zero norm 2nd order AC matrix on the fly, and generate the
                  % "out" condition matrix on the fly to conserve memory
                  
                  % decorr wrt specific attentional stuff
                  B=uu'*SAatt(:,:,:,resampidx,inidx,1)*uu;
                  eigHatt(:,:,resampidx,inidx)=...
                      B^-1 * eigSatt(:,:,resampidx,inidx);
                  
                  B=uu'*SAatt(:,:,:,resampidx,inidx,2)*uu;
                  eigHoutatt(:,:,resampidx,inidx)=...
                      B^-1 * eigSoutatt(:,:,resampidx,inidx);
               end
            end      % for resampidx
         end     % for inidx
         
         fprintf('T2: ');
         eigrange=10.^(linspace(7,0,ktestcount));
         
         %keyboard
         
         for rr=1:resptouse,
            for att1=1:attdocount,
               for eigidx=1:ktestcount,
                  
                  %tsigidx=find(sigidx(1:eigidx,rr+respbase));
                  tsigidx=1:eigidx;
                  elen=length(tsigidx);
                  
                  %compute T2 of difference between eigH and eigHout
                  % X is PC x resamp x in/out cond
                  
                  if ~DOREGNORM & elen>0,
                     % do sharp cutoff T2 test
                     
                     % pull out appropriate pair of kernels to compare
                     if DOPAIRS,
                        a1=attpairs(att1,1);
                        a2=attpairs(att1,2);
                        X=cat(4,eigHatt(tsigidx,:,:,a1),...
                              eigHatt(tsigidx,:,:,a2));
                     else
                        X=cat(4,eigHatt(tsigidx,rr+respbase,:,att1),...
                              eigHoutatt(tsigidx,rr+respbase,:,att1));
                     end
                     
                     % isolate latency range of interest (for time
                     % kernels) or time bin of interest (for PFTH
                     % kernels)
                     if respfmtcode==0,
                        tmax=min([maxlag(2) 7]);
                        trange=-maxlag(1)+(1:tmax);
                        tlen=length(trange);
                        X=reshape(X(:,trange,:,:),elen*tlen,resampcount,2);
                     else
                        X=reshape(X(:,rr+respbase,:,:),elen,resampcount,2);
                     end
                     
                     % Hoteling's T algorithm from Good. Pertmutation
                     % Tests. p. 65
                     E=mean(X(:,:,1),2)-mean(X(:,:,2),2);
                     mX=mean(X,2);
                     X=X-repmat(mX,[1 resampcount 1]);
                     V=(X(:,:,1)*X(:,:,1)' + X(:,:,2)*X(:,:,2)') ./ ...
                       (resampcount-2);
                     
                     % if V were nice
                     Vinv=V^-1;
                     
                     % pseudo inverse... seems to flail in extremes
                     %Vinv=svdinv(V,0.000001);
                     
                     % ignore correlations between channels?
                     % amplifies noise but may be the only stable option
                     %sV=1./diag(V);
                     %sV(find(sV>median(sV)*10^6))=0;
                     %Vinv=diag(sV);
                     
                     T2(att1,sampidx,eigidx,rr)=E' * Vinv * E ./ ...
                         size(X,1);
                  elseif elen>0,
                     
                     % this is old and non-useful
                     sadj=1./(1+eigrange(eigidx)./sd(:).^2);
                     sadj=repmat(sadj,[1 resampcount]);
                     
                     ein=squeeze(eigHatt(:,rr+respbase,:,att1));
                     eout=squeeze(eigHoutatt(:,rr+respbase,:,att1));
                     
                     X=cat(3,uu*(ein.*sadj),uu*(eout.*sadj));
                     
                     X=X./max(X(:));
                     E=mean(X(:,:,1),2)-mean(X(:,:,2),2);
                     mX=mean(X,2);
                     X=X-repmat(mX,[1 resampcount 1]);
                     V=(X(:,:,1)*X(:,:,1)' + ...
                        X(:,:,2)*X(:,:,2)')./(resampcount-2);
                     %[u,s,v]=svd(V);
                     %s(find(s<s(1,1)/10^5))=0;
                     %s(find(s>0))=1./s(find(s>0));
                     %Vinv=v*s*u';
                     %Vinv=V^-1;
                     sV=1./diag(V);
                     sV(find(sV>median(sV)*10^6))=0;
                     Vinv=diag(sV);
                     T2(sbase+att1,eigidx,rr)=E' * Vinv * E ./spacecount;
                     
                  end
                  
               end %for eigidx
            end
            fprintf('.');
         end
         fprintf('\n');
         
         % save average across resamples
         eigS(:,:,sbase+(1:(attuse-1)))=squeeze(mean(eigSatt,3));
         eigSout(:,:,sbase+(1:(attuse-1)))=squeeze(mean(eigSoutatt,3));
         eigH(:,:,sbase+(1:(attuse-1)))=squeeze(mean(eigHatt,3));
         eigHout(:,:,sbase+(1:(attuse-1)))=squeeze(mean(eigHoutatt,3));
         
      end % if inidx==attuse-1
      
   end  % for attidx
   
   if PRELOADMOV==2,
      delete(tstimfile);
   end
end

for attidx=1:attdocount,
   for eigidx=1:ktestcount,
      for rr=1:resptouse,
         tt=T2(:,2:end,eigidx,rr);
         tm=T2(attidx,1,eigidx,rr);
         tt=sort([tt(:); tm]);
         if ~isnan(tm),
            p(attidx,eigidx,rr)=...
                1-min(find(tm<=tt)-1)/(length(tt));
         else
            disp('NAN TM!!!!!!!!!!!!!!');
         end
      end
   end
end

if 0,
% ok, we've loaded eigH for all resample conditions for
% this attentional state. now do the T2 test

disp('hoteling''s T time...');
resptouse=4;  % don't use all respcount eigH's
respbase=4;
latidx=6;
T2=zeros(attcount,spacelim,resptouse);
p=ones(attuse,spacelim,resptouse);

for respidx=1:resptouse,
   for attidx=2:attcount,
      for eigidx=1:spacelim,
         
         %compute T2 of difference between eigH and eigHout
         % X is PC x resamp x in/out cond
         X=reshape(cat(4,eigH(1:eigidx,respidx+respbase,:,attidx),...
                       eigHout(1:eigidx,respidx+respbase,:,attidx)),...
                   eigidx,resampcount,2);
         E=mean(X(:,:,1),2)-mean(X(:,:,2),2);
         mX=mean(X,2);
         X=X-repmat(mX,[1 resampcount 1]);
         V=(X(:,:,1)*X(:,:,1)' + X(:,:,2)*X(:,:,2)')./(resampcount-2);
         Vinv=V^-1;
         %if 1| eigidx<resampcount,
         %   Vinv=V^-1;
         %else
         %   [u,s,v]=svd(V);
         %   sinv=1./diag(s);
         %   sinv(resampcount:end)=0;
         %   Vinv=v'*diag(sinv)*u;
         %end
         T2(attidx,eigidx,respidx)=E' * Vinv * E;
      end %for eigidx
      
      if mod(attidx,15)==0,
         fprintf('.');
      end
      if exist('BATQUEUEID') & BATQUEUEID>0,
         dbsetqueue(BATQUEUEID,1000+attidx);
      end
   end
end
fprintf('\n');

for attidx=2:attuse,
   for eigidx=1:spacelim,
      for respidx=1:resptouse,
         tt=sort(T2([attidx (attuse+1):end],eigidx,respidx));
         p(attidx,eigidx,respidx)=...
             1-min(find(T2(attidx,eigidx,respidx)<=tt)-1)/(noisecount+1);
      end
   end
end

eigH=squeeze(mean(eigH,3));
eigHout=squeeze(mean(eigHout,3));
end

meigH=squeeze(mean(eigH(:,:,(attuse+1):end),3));
seigH=squeeze(std(eigH(:,:,(attuse+1):end),1,3)) .* ...
      sqrt((noisecount-1)/noisecount);
meigHout=squeeze(mean(eigHout(:,:,(attuse+1):end),3));
seigHout=squeeze(std(eigHout(:,:,(attuse+1):end),1,3)) .* ...
         sqrt((noisecount-1)/noisecount);
seigHdiff=squeeze(std(eigH(:,:,(attuse+1):end)-...
                      eigHout(:,:,(attuse+1):end),1,3)) .* ...
          sqrt((noisecount-1)/noisecount);

%
% FIGURE OUT WHAT THIS HAS TO TO WITH THE TARGETS!
%

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
   
   % transform patches into PC space
   %tu=diag(1./sd(:,1)) * uu';
   teigpatches=uu' * mpatches;
   beigpatches=uu' * stim';
   eigpatches=uu' * fpatches;
   
   % compute mean and std across patches in PC space
   meigpatches=squeeze(mean(eigpatches,2));
   mbeigpatches=squeeze(mean(beigpatches,2));
   seigpatches=squeeze(std(eigpatches,1,2));   
   sbeigpatches=squeeze(std(beigpatches,1,2));   
   
   attcc=zeros(spacelim,attuse,3);
   ut=zeros(spacecount,spacelim,attuse,3);
   for attidx=1:attuse,
      % scale each stimulus eigenvector by kernel coefficient
      ut(:,:,attidx,1) = uu * diag(eigH(:,latidx,attidx));
      ut(:,:,attidx,2) = uu * diag(eigHout(:,latidx,attidx));
   end
   
   % take cumulative sum to get kernel at successively greater
   % eigenvector cutoffs
   uts=cumsum(ut,2);
   
   targcount=size(mpatches,2);
   stdpatches=std(fpatches,1,2);
   stdpatches(find(stdpatches==0))=1;
   tpatches=mpatches-repmat(mean(fpatches,2),[1 targcount]);
   tpatches=tpatches./repmat(stdpatches,[1 targcount]);
   
   for eigidx=1:spacelim,
      for attidx=1:attuse,
         % corr of attin kernel with target
         if sum(abs(uts(:,eigidx,attidx,1))) & sum(abs(tpatches(:,attidx))),
            xc=corrcoef(tpatches(:,attidx),uts(:,eigidx,attidx,1));
            attcc(eigidx,attidx,1)=xc(2,1);
         end
         % corr of attout kernel with target
         if sum(abs(uts(:,eigidx,attidx,2))) & sum(abs(tpatches(:,attidx))),
            xc=corrcoef(tpatches(:,attidx),uts(:,eigidx,attidx,2));
            attcc(eigidx,attidx,2)=xc(2,1);
         end
         % corr of mean kernel with target
         if sum(abs(uts(:,eigidx,1,1))) & sum(abs(tpatches(:,attidx))),
            xc=corrcoef(tpatches(:,attidx),uts(:,eigidx,1,1));
            attcc(eigidx,attidx,3)=xc(2,1);
         end
      end
   end
   
   % measure similarity of different targets
   if DOPAIRS,
      % for each target pair, incrementig included eigs, and norm
      % or not normalized
      targsim=zeros(attdocount,spacelim,2);
      for attidx=1:attdocount,
         for eigidx=1:spacelim,
            targsim(attidx,eigidx,1)= ...
                sum(tpatches(1:eigidx,attpairs(attidx,1)+1) .* ...
                    tpatches(1:eigidx,attpairs(attidx,2)+1));
            targsim(attidx,eigidx,2)= targsim(attidx,eigidx,1) ./ ...
                (sqrt(sum(tpatches(1:eigidx,attpairs(attidx,1)+1).^2)).*...
                 sqrt(sum(tpatches(1:eigidx,attpairs(attidx,2)+1).^2)));
         end
      end
   end
else
   % do something to figure out target numbers here!
   targsim=zeros(attdocount,spacelim,2);
   attcc=zeros(spacelim,attuse,3);
   
end

% clear various input and intermediate stuff that wastes space
clear stim resp bstim bresp z data
clear tSR tsSA1 tsSA2 ttSA
clear sSA1 sSA2 tSA
clear *att SA SR V Vinv X B timage xc u s v tn tt tu

if spacecount>10,
   % clear unwanted space wasters
   clear H ttSA tsSA1 tsSA2 tSR
   clear tSA sSA0 sSA1 sSA2 sSA2out SR
end

return
 


         tgoodidx=find(~isnan(resp(:,respidx)));
         plotrange=1:min([200 length(tgoodidx)]);
         for xx=1:spacecount,
            subplot(attuse,spacecount,(attidx-1)*spacecount+xx);
            scatter(stim(tgoodidx(plotrange),xx),...
                    resp(tgoodidx(plotrange),respidx));
            title(sprintf('stim %d att %d',xx,attidx));
            
            [p,s]=polyfit(stim(tgoodidx,xx),resp(tgoodidx),1);
            x=unique(stim(tgoodidx));
            hold on
            plot(x,p(1).*x + p(2),'k-');
            hold off
         end

         





tbincount=respcount;
stylestr={'k-','c--','b--','r--','g--'};

figure(1);
clf
for sfsidx=1:sfscount,
   
   mm=max(max(max(abs(mH(:,:,sfsidx,1,:)))));
   
   for tidx=1:tbincount,
      subplot(sfscount,tbincount,tidx+(sfsidx-1)*tbincount);
      hold on
      for attidx=1:attcount,
         errorbar(1:spacecount,mH(:,tidx,sfsidx,1,attidx),...
                  eH(:,tidx,sfsidx,1,attidx)*2,stylestr{attidx});
      end
      hold off
      title(sprintf('lat=%d sfs=%d',tidx,sfsidx));
      axis([1 spacecount -mm mm]);
      
   end
end
legend('A','1','2','3','4');

