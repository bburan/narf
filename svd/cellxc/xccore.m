% xccore.m
%
% script to compute unbiased kernel from preprocessed stimulus and
% response matrices (ie, stim T x X and resp T x M supplied)
% 
% returns:
%
% mean and stderr on kernel estimates: (X x tbincount x M x sfscount)
% mH=mean(H,5);
% eH=std(H,1,5) .* sqrt((resampcount-1)/resampcount);
% 
% mean and stderr on raw (biased) XC: (X x tbincount x M x sfscount)
% mSR=mean(SR,4);
% eSR=std(SR,1,4).* sqrt((resampcount-1)/resampcount);
% 
% % mean stim/response across resamples
% mSall=mean(mS,3);
% mRall=mean(mR,2);
% ntot=n;
%
% created SVD 1/23/03 - ripped off of cellxc.m
%

disp('xccore.m:');

%
% define resampling regimes here
%
% option to branch according to resampfmt
if params.resampcount==1 | length(resp)==0,
   % bootstrap resampling
   rstartidx=1;
   rendidx=size(resp,1);
elseif params.resampfmt==1,
   [rstartidx,rendidx]=resampsegs(resp,params.resampcount);
   
elseif params.resampfmt==2,
   % shuffled response resampling
   rstartidx=starttimes(fidx,attidx);
   rendidx=stoptimes(fidx,attidx);
   resp=resampnoise(resp,params.resampcount);
   rsize=[rsize,size(resp,3)];
   resp=resp(:,:); % reshape to 2D matrix
   
elseif params.resampfmt==3,
   % gamma/poisson spiking model resampling   
   rstartidx=1;
   rendidx=size(resp,1);
   respcount=1;
   rsize(2)=1;
   [mu,alpha,beta]=reversepoisson(resp);
end

%keyboard
%stim(:,spacecount+1)=1;
%spacecount=spacecount+1;

if firstseg,
   disp('initializing matrices');
   if 0 & (diff(params.maxlag)==0 | rsize(2)==1),
      corrmtxcount=1;
      singlesSA=1;
   else
      corrmtxcount=rsize(2);
      singlesSA=0;
   end
   
   % first time, set kernel & ac matrices to zero
   SR=zeros(spacecount,diff(params.maxlag)+1,respcount,params.resampcount);
   n=zeros(respcount,params.resampcount);
   mS=zeros(spacecount,respcount,params.resampcount);
   mR=zeros(respcount,params.resampcount);
   tSA=zeros(diff(params.maxlag)*2+1,corrmtxcount,params.resampcount);
   sSA2=zeros(spacecount,spacecount,corrmtxcount,params.resampcount);
   if params.decorrspace==4,
      sSAfull=zeros(spacecount*(diff(params.maxlag)+1),spacecount*(diff(params.maxlag)+1),...
          params.resampcount);
   end
   
   % global stimulus temporal autocorr. works better than resamped?
   tSA1=zeros(diff(params.maxlag)*2+1,corrmtxcount,params.resampcount);
   nSA1=0;
   firstseg=0;
end


fprintf('resamp=');
for resampidx=1:params.resampcount,
   fprintf('%d ',resampidx);
   for respidx=1:respcount,
      if params.resampfmt~=3,
         tr=resp(:,respidx);
         tr([1:rstartidx(resampidx,respidx)-1 ...
             rendidx(resampidx,respidx)+1:end],:)=nan;
         %tr(rstartidx(resampidx,respidx):rendidx(resampidx,respidx),:)=nan;
      else
         tr=reversepoisson(resp(:,2:end),alpha,beta);
      end
      
      rgoodidx=find(~isnan(tr));
      tmR=sum(resp(rgoodidx,respidx))';
      tmS=sum(stim(rgoodidx,:))';
      tn=length(rgoodidx);
      
      tSR=zeros(spacecount,diff(params.maxlag)+1);
      
      for tt=params.maxlag(1):params.maxlag(2),
         %fprintf('.');
         trg=rgoodidx;
         trg=trg(find(trg-tt>0 & trg-tt<size(stim,1)));
         tSR(:,tt-params.maxlag(1)+1)=(resp(trg,respidx)'*stim(trg-tt,:))';
      end
      
      if respidx<=corrmtxcount & length(rgoodidx)>0,
         if params.decorrspace==2,
            tsSA2=stim(rgoodidx,:)'*stim(rgoodidx,:);
            %fprintf('*');
         elseif params.decorrspace==4 & respidx==1,
            tbincount=diff(params.maxlag)+1;
            stimplus=zeros(length(rgoodidx),spacecount*tbincount);
            for tbidx=1:tbincount,
                stimplus(:,(1:spacecount)+(tbidx-1)*spacecount)=...
                    stim(rgoodidx+tbidx-params.maxlag(2)-1,:);
            end
            tsSAfull=stimplus'*stimplus;
            clear stimplus
         end
         ttSA=zeros(diff(params.maxlag)*2+1,1);
         for xx=1:spacecount,
            tstim=stim(rgoodidx,xx);
            ttSA=ttSA+xcorr(tstim,diff(params.maxlag),'biased') ./ ...
                 spacecount.*length(rgoodidx);
         end
         %fprintf('*\n');
      end
      
      for rr=1:params.resampcount,
         % add means to all resamp channels. this is to
         % make it so that everything is centered around
         % the same DC!
         if params.resampcount==1 | rr~=resampidx,
            % otherwise add outputs to running total
            SR(:,:,respidx,rr)=SR(:,:,respidx,rr)+tSR;
            mS(:,respidx,rr)=mS(:,respidx,rr)+tmS;
            mR(respidx,rr)=mR(respidx,rr)+tmR;
            n(respidx,rr)=n(respidx,rr)+tn;
            tSA(:,respidx,rr)=tSA(:,respidx,rr)+ttSA;
            
            if params.decorrspace==2 & respidx<=corrmtxcount,
               sSA2(:,:,respidx,rr)=sSA2(:,:,respidx,rr)+tsSA2;
            elseif params.decorrspace==4 & respidx==1,
               sSAfull(:,:,rr)=sSAfull(:,:,rr)+tsSAfull;
            end
         end
      end
   end
   
   % update queue if active
   %dbsetqueue;
end
fprintf('\n');

if diff(params.maxlag)>1,
   for xx=1:spacecount,
      
      tstim=stim(:,xx)-mean(stim(:,xx));
      tSA1=tSA1+repmat(xcorr(tstim,diff(params.maxlag),'biased')./ ...
                       spacecount.*size(stim,1),...
                       [1 corrmtxcount params.resampcount]);
   end
else
   tSA1(:)=tSA1(:) + size(stim,1);
end
nSA1=nSA1+size(stim,1);

% zeroth order normalization. ie, divide sums by the number of
% samples to get appropriate means
[SR,mS,mR,tSA,sSA2]=normalize0(SR,n,mS,mR,tSA,sSA2,params.meansub);
tSA1=tSA1./nSA1;

fprintf('params.decorrspace=%d\n',params.decorrspace);

% special case preparation for decorrelation. generally use
% decorrspace=2, which doesn't require anything
if params.decorrspace==3,
   load([params.stimfiles{fidx},'.AC.mat'],'sSA2');
elseif params.decorrspace==1,
   
   % normalize by power only... ie, force corr matrix to be
   % diagonal. dunno when this would be useful
   for ii=1:size(sSA2,3).*size(sSA2,4),
      sSA2(:,:,ii)=diag(diag(sSA2(:,:,ii)));
   end
   
elseif 0,
   disp('taking SINGLE sSA2!');
   % take average sSA2 to speed things up and possibly reduce
   % noise in kernel estimates... hell, it worked in V4
   sSA2=mean(sSA2,4);
end

% remove off peak terms of temporal correlation if no temporal
% decorr requested
if isfield(params,'decorrtime') & ~params.decorrtime,
   disp('skipping temporal normalization');
   tSA1(1:diff(params.maxlag),:,:)=0;
   tSA1((1:diff(params.maxlag))+diff(params.maxlag)+1,:,:)=0;
end

mSA2=mean(sSA2(:,:,:),3);
powbiased=mSA2*ones(spacecount,1);

% 8/7/02: use single tSA for all resamples
if spacecount==1,
   
   % small spatial dimension case means shouldn't do spatial decorr
   [H,lambda,pctvar]=normalizetime(SR,tSA,params.sfscount);
   
elseif params.decorrspace==4,
   keyboard
   [H,lambda,pctvar]=normalizefull(SR,sSAfull,params.sfscount,...
                                   params.sfsstep,1,params.smoothtime);
   
elseif params.sharpspacenorm==1,
   
   % standard sharp edged pseudo-inverse
   [H,lambda]=normalize(SR,sSA2,tSA1,...
                        10.^(linspace(-1,-params.sfsstep,params.sfscount)));
   powunbiased=normalize(powbiased,mSA2,[],...
                         10.^(linspace(-1,-params.sfsstep,params.sfscount)));
   
elseif params.sharpspacenorm==2,
   % smooth pseudo-inverse cutoff
   [H,lambda,pctvar]=normalizeshr(SR,sSA2,tSA1,params.sfscount,...
                     params.sfsstep,1,params.smoothtime); % topSR=1
   powunbiased=normalizeshr(powbiased,mSA2,[],params.sfscount,...
                            params.sfsstep,1,0);
else
   % smooth pseudo-inverse cutoff
   [H,lambda,pctvar]=normalizereg(SR,sSA2,tSA1,params.sfscount,...
                     params.sfsstep,1,params.smoothtime); % topSR=1
   powunbiased=normalizereg(powbiased,mSA2,[],params.sfscount,...
                            params.sfsstep,1);
   if isfield(params,'fixresbias') & params.fixresbias,
      for sfsidx=1:params.sfscount,
         fff=find(abs(powunbiased(:,sfsidx))>1);
         
         pp=1-(1-abs(powunbiased(:,sfsidx))).^2;
         pp(abs(powunbiased(:,sfsidx))>1)=1;
         fff=1:length(pp);
         
         if length(fff)>0,
            H(fff,:,sfsidx,:,:)=H(fff,:,sfsidx,:,:) ./ ...
                repmat(abs(powunbiased(fff,sfsidx)./pp),...
                       [1 diff(params.maxlag)+1 1 ...
                        respcount params.resampcount]);
            powunbiased(fff,sfsidx)=pp.*sign(powunbiased(fff,sfsidx));
         end
      end
   end
   
   if params.sffiltsmooth,
      %if size(SR,2)>1 & params.decorrtime,
      %   [tH,lambda]=normalizesmooth(SR,tSA1,0.001);
      %   tH=permute(tH,[1 2 4 5 6 3]);
      %else
      %   tH=SR;
      %end
      %
      %[H,lambda]=normalizesmooth(tH,sSA2,params.sfscount,params.sfsstep);
      %powunbiased=normalizesmooth(powbiased,mSA2,params.sfscount,...
      %                            params.sfsstep);
      
      disp('smoothing or=10/sf=6...');
      %keyboard
      H=kernsmoothorsf2(H,10,6,params.kernfmt);
   end
end

irS=squeeze(std(SR,1,1));
irH=squeeze(std(H,1,1));

allt=std(reshape(permute(H,[1 2 4 5 3]),spacecount*(diff(params.maxlag)+1)* ...
                 respcount*params.resampcount,params.sfscount));
alls=std(reshape(permute(H,[1 3 4 5 2]),spacecount*params.sfscount* ...
                 respcount*params.resampcount,(diff(params.maxlag)+1)));
tH=squeeze(std(reshape(permute(H,[1 4 5 2 3]),spacecount* ...
                       respcount*params.resampcount,...
                       (diff(params.maxlag)+1),params.sfscount)));

mH=mean(H,5);
eH=std(H,1,5) .* sqrt(params.resampcount-1);

mSR=mean(SR,4);
eSR=std(SR,1,4).* sqrt(params.resampcount-1);

% take mean stim/response across resamples
mSall=mean(mS,3);
mRall=mean(mR,2);
ntot=n;









