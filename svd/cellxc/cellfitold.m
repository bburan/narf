% cellfit.m
%
% now it's time to fit the data!!!  ack!   ????

% we have two things to fit: the svd cutoff (sfsidx) and the
% shrinkage filter fudge factor (sfiltsigma)

DOTHRESH=0;
spacecount=size(mH,1);
if diff(maxlag)>0,
   tbincount=min([maxlag(2)+1 14]);
else
   tbincount=size(mH,2);
end
sfscount=size(mH,3);
respcount=size(mH,4);
attcount=size(mH,5);
lagrange=(1-maxlag(1)):(tbincount-maxlag(2)+10);

sfsidxtouse=(1:sfscount)';  % test all sfs cutoffs

% define a range of sigmas for combining resampled estimates (via
% shrinkage filter)
%sigrange=[1.0 1.2 1.4 1.7 2.0 2.4 2.8 3.3 3.8 4.6 5.2 6.0 7.0 8.0
%9.0];
if respfmtcode==0,
   if sffiltsigma==1,
      sigrange=4;
   else
      %sigrange=exp(linspace(1,2.5,sffiltsigma));
      sigrange=exp(linspace(log(3),log(10),sffiltsigma));
      %sigrange=3.^(linspace(1,2.1,sffiltsigma));
      %sigrange=(linspace(sqrt(3),sqrt(12),sffiltsigma)).^2;
   end
   
else
   sigrange=exp(linspace(0,1.5,sffiltsigma));
   %sigrange=exp(linspace(-1,log(2.5),10));
end
if 0, % run fast  for debugging
   sigrange=exp([0 1]);
   sfsidxtouse=0:2:10;
   sfscount=length(sfsidxtouse);
end
sigcount=length(sigrange);

if ~exist('nloutparm'),
   nloutparm=1;
end
nlcount=3;

% for PFTH preds, separate score for each latency bin
if respfmtcode==0,
   latcount=1;
else
   latcount=tbincount;
end
xc=zeros(sfscount,sigcount,latcount,nlcount,attcount);
mutinfo=zeros(sfscount,sigcount,latcount,nlcount,attcount);
hingeparms=zeros(3,latcount,nlcount,attcount);
threshparm=zeros(sfscount,sigcount,latcount,nlcount,attcount);
tgoodlen=zeros(latcount,attcount);
sfsfit=zeros(latcount,nlcount,attcount);
sigfit=zeros(latcount,nlcount,attcount);
threshfit=zeros(nloutparm,latcount,nlcount,attcount);
xcnl=zeros(nloutparm,latcount,nlcount,attcount);
hf=zeros(spacecount,tbincount,nlcount,attcount);

%samplecount=6; %number of evenly distributed sample kernels to save
%hfsample=zeros(spacecount,tbincount,samplecount,attcount);

for attidx=1:attcount,
   fprintf('Fitting for attentional state %d.\n',attidx);
   
   % load stim and response segments to fit. although this will
   % generally be the same file for all attentional states, it
   % theoretically could be different. thus reload the stimulus and
   % response for each attentional state
   
   % adjust fitstart frame with extra frames to beginning so that
   % first valid prediction bin lines up with the beginning of fit
   % data.
   fitfilecount=sum(fitfile(:,attidx)>0);
   stim=[];
   resp=[];
   for fitidx=1:fitfilecount,
      if respfmtcode==0,
         leadbincount=min([fitstartframe(fitidx,attidx)-1 tbincount-1]);
      else
         leadbincount=0;
      end
      tfitstartframe=fitstartframe(fitidx,attidx)-leadbincount;
      
      % if scale_pix not defined, use crf parms
      tstimloadparms=stimloadparms;
      if strcmp(stimloadcmd,'loadimfile') & ...
            length(stimloadparms)>0 & stimloadparms{1}==0,
         tstimloadparms{1}=stimloadparms{3} .* ...
             stimcrfs(fitfile(fitidx,attidx))./stimwindowcrf;
      end
      tstim=feval(stimloadcmd,stimfiles{fitfile(fitidx,attidx)},...
                  tfitstartframe,fitstopframe(fitidx,attidx),...
                  tstimloadparms{:});
      
      % filter stimulus segment if selected
      if ~isempty(stimfiltercmd),
         tstim=feval(stimfiltercmd,tstim,stimfilterparms{:});
      end
      iconside=size(tstim);
      tmovlen=iconside(end);
      iconside=iconside(1:(end-1));
      if length(size(tstim))>2,
         % reshape all spatial dims into one
         tstim=reshape(tstim,prod(iconside),tmovlen);
      end
      
      stim=cat(2,stim,tstim);
      
      tresp=feval(resploadcmd,respfiles{fitfile(fitidx,attidx)},...
                  resploadparms{:});
      tresp=tresp(tfitstartframe:fitstopframe(fitidx,attidx),:,:);
      
      % filter response (eg smooth, pick out ICs, PCs? currently not used)
      if ~isempty(respfiltercmd),
         tresp=feval(respfiltercmd,tresp,respfilterparms{:});
      end
      predlen=fitstopframe(fitidx,attidx)-tfitstartframe+1;
      
      % mark first bins invalid so they don't get used in pred eval
      tresp(1:tbincount,:,:)=nan;
      
      resp=cat(1,resp,tresp);
   end
   
   % figure out actual valid range of response being predicted
   tr=resp(:,:,attidx);
   
   % smooth actual psth if movformat.predsmoothsigma > 0
   if batchdata.predsmoothsigma > 0 & respfmtcode==0,
      %& cellfiledata(fitfile).repcount==1,
      %fprintf('Smoothing actual response with optimal filter...\n');
      %pfilt=[1/9 2/9 1/3 2/9 1/9]';
      
      fprintf('Smoothing fit response sigma=%.1f bins...\n',...
              batchdata.predsmoothsigma);
      tt=(-10:10)';
      pfilt=exp(-tt.^2/(2*batchdata.predsmoothsigma.^2))./...
            (sqrt(2*pi)*batchdata.predsmoothsigma);
      pfilt=pfilt(find(pfilt > 0.1*max(pfilt)));
      fr=tr;
      fr(isnan(fr))=0;
      fr=conv2(fr,pfilt,'same');
      tr(~isnan(tr))=fr(~isnan(tr));
   end
   
   % normalize valid bins of actual response
   tgoodix={};
   for r=1:latcount,
      tgoodidx{r}=find(~isnan(tr(:,r)));
      tgoodlen(r,attidx)=length(tgoodidx{r});
   end
   
   % loop through each attentional state, sfs cutoff and sigma fudge
   % factor and predict the response to the fit data
   for r=1:respcount,
      if meansub,
         % subtract exploratory mean from stimulus to get DC level
         % set correctly
         tstim=stim-repmat(mSall(:,r,attidx),[1 size(stim,2)]);
      else
         tstim=stim;
      end
      for sfsidx=1:sfscount,
         % smooth the kernel??!?!?!
         smm=mH(:,(1:tbincount)-maxlag(1),sfsidx,r,attidx);
         if resampcount>1,
            sms=eH(:,(1:tbincount)-maxlag(1),sfsidx,r,attidx);
         end
         
         for sigidx=1:sigcount,
            fprintf('Testing fit (att=%d sfsidx=%d sig=%.2f r=%d)...\n',...
                    attidx,sfsidxtouse(sfsidx),sigrange(sigidx),r);
            % current kernel is shrinkage-filtered from current (sfs,r,attidx)
            % only use causal time lags for pred!
            if resampcount>1,
               smd=abs(smm)./(sms+(sms==0)) ./ sigrange(sigidx);
               if ~DOTHRESH, % old--shrinkage filter
                  smd=(1-smd.^(-2));
                  smd=smd.*(smd>0);
                  smd(find(isnan(smd)))=0;
                  thf=smm.*smd;
               else
                  % new -- threshold by # of std errs
                  thf=smm.*(smd>1);
               end
            else
               thf=smm;
            end
            
            for nlidx=1:nlcount,
               if respfmtcode==1,
                  
                  % don't sum over latency dimension... this returns a
                  % lat X 1 X latency matrix
                  tmod_psth=kernpredict(thf,tstim,1,0,0);
                  %tmod_psth=kernpredict(thf,tstim,spacecount,1,0);
                  tmod_psth=squeeze(tmod_psth);
                  
                  % "resp" space comes back in columns. time still in
                  % rows. need to adjust normalization routines to
                  % deal with more than one column in tmod_psth
               elseif nlidx==1,    % no rectification
                  seplinpred=kernpredict(thf,tstim,spacecount,0);
                  linpred=sum(seplinpred,2);
                  tmod_psth=linpred;
               elseif nlidx==2,    % rectify summed output
                  %tmod_psth=kernpredict(thf,tstim,1,1);
                  % to speed things up, don't re-predict. just
                  % rectify previously generated linear prediction
                  tmod_psth=linpred.*(linpred>0);
               elseif nlidx==3,    % rectify each spatial channel individually
                  %tmod_psth=kernpredict(thf,tstim,spacecount,1);
                  tmod_psth=seplinpred;
                  tmod_psth(find(tmod_psth<0))=0;
                  tmod_psth=sum(tmod_psth,2);
               elseif nlidx==4,
                  tmod_psth=linpred;
                  [t0,res]=findthresh(tmod_psth(tgoodidx{1}),tr(tgoodidx{1}));
                  threshparm(sfsidx,sigidx,1,nlidx,attidx)=t0;
                  tmod_psth(tgoodidx{1})=thresh(t0,tmod_psth(tgoodidx{1}));
               elseif nlidx==5,
                  tmod_psth=seplinpred;
                  [t0,res]=findthresh(tmod_psth(tgoodidx{1},:),...
                                      tr(tgoodidx{1}));
                  threshparm(sfsidx,sigidx,1,nlidx,attidx)=t0;
                  tmod_psth=thresh(t0,tmod_psth);
               end
               
               % remove leading bins that don't matter
               for ridx=1:latcount,
                  % normalize for appropriate correlation
                  % coefficient calculation
                  if std(tmod_psth(tgoodidx{ridx},ridx))>0,
                     xc(sfsidx,sigidx,ridx,nlidx,attidx)=...
                         xcov(tr(tgoodidx{ridx},ridx),...
                              tmod_psth(tgoodidx{ridx},ridx),0,'coeff');
                     %coh=cohere(tr(tgoodidx{ridx},ridx),...
                     %           tmod_psth(tgoodidx{ridx},ridx),...
                     %           32,70,[],16);
                     %mutinfo(sfsidx,sigidx,ridx,nlidx,attidx)=...
                     %    -sum(log2(1-coh.^2));
                  else
                     disp('zero pred cond');
                  end
               end
            end
         end % for sigidx

         if not(isempty(BATQUEUEID)),
            % record latest (sfsidx,attidx) predicted
            dbsetqueue(BATQUEUEID,sfsidx+(attidx-1)*100+1000);
         end
      end     % for sfsidx
   end    % for r .. respidx
  
   % generate optimal kernels and save in hf
   for nlidx=1:nlcount,
      for latidx=1:latcount,
         
         % figure out max xc
         txc=xc(:,:,latidx,nlidx,attidx);
         
         %smooth txc?
         if size(txc,1)>=8 & size(txc,2)>=3,
            [XX,YY]=meshgrid(-1:1,-1:1);
            sx=0.5; sy=1.0;
            % fvvs
            % sx=1.5; sy=2.0;
            g=exp(-XX.^2./(2*sx.^2) - YY.^2./(2*sy.^2));
            g=g./sum(g(:));
            txc=rconv2(txc,g);
         end
         maxidx=min(find(txc==max(txc(:))));
         [sfsfit(latidx,nlidx,attidx),sigfit(latidx,nlidx,attidx)]=...
             ind2sub([sfscount,sigcount],maxidx);
         if nlidx>=4,
            threshfit(latidx,nlidx,attidx)=...
                threshparm(sfsfit(latidx,nlidx,attidx),...
                           sigfit(latidx,nlidx,attidx),latidx,nlidx,attidx);
         end
         fprintf('attidx=%d nlidx=%d latidx=%d: (sfs,sig)=(%d,%d)\n',...
                 attidx,nlidx,latidx,sfsfit(latidx,nlidx,attidx),...
                 sigfit(latidx,nlidx,attidx));
      end
      
      % fill hf with best fit kernel according to sigfit and sfsfit
      %
      % use same fit values for all latencies (in respfmtcode==1)
      %
      sfsfituse=floor(median(sfsfit(:,nlidx,attidx)));
      sigfituse=floor(median(sigfit(:,nlidx,attidx)));
      
      for latidx=1:latcount,
         if latcount>1,
            smm=mH(:,latidx,sfsfituse,1,attidx);
            if resampcount>1,
               sms=eH(:,latidx,sfsfituse,1,attidx) .* sigrange(sigfituse);
               smd=abs(smm)./(sms+(sms==0));
               if ~DOTHRESH,
                  smd=1-smd.^(-2);
                  smd=smd.*(smd>0);
                  smd(find(isnan(smd)))=0;
                  hf(:,latidx,nlidx,attidx)=smm.*smd;
               else
                  hf(:,latidx,nlidx,attidx)=smm.*(smd>1);
               end
            else
               hf(:,latidx,nlidx,attidx)=smm;
            end
         else
            smm=mH(:,(1:tbincount)-maxlag(1),sfsfituse,1,attidx);
            if resampcount>1,
               sms=eH(:,(1:tbincount)-maxlag(1),sfsfituse,1,attidx) .* ...
                   sigrange(sigfituse);
               smd=abs(smm)./(sms+(sms==0));
               if ~DOTHRESH,
                  smd=1-smd.^(-2);
                  smd=smd.*(smd>0);
                  smd(find(isnan(smd)))=0;
                  hf(:,:,nlidx,attidx)=smm.*smd;
               else
                  hf(:,:,nlidx,attidx)=smm.*(smd>1);
               end
            else
               hf(:,:,nlidx,attidx)=smm;
            end
         end
      end
      
      % fit output nonlinearities for each hf
      for nloutidx=1:nloutparm,
         if nloutidx==1,
            % no rectification
            seplinpred=kernpredict(hf(:,:,nlidx,attidx),tstim,spacecount,0);
            linpred=sum(seplinpred,2);
            tmod_psth=linpred;
            t0=min(tmod_psth(:));
         elseif nloutidx==2,
            % rectify summed output
            tmod_psth=linpred;
            [t0,res]=findthresh(tmod_psth(tgoodidx{1}),tr(tgoodidx{1}));
         elseif nloutidx==3,
            % rectify each spatial channel
            tmod_psth=seplinpred;
            [t0,res]=findthresh(tmod_psth(tgoodidx{1},:),...
                                tr(tgoodidx{1}));
         end
         
         tmod_psth=thresh(t0,tmod_psth);
         threshfit(nloutidx,1,nlidx,attidx)=t0;
         xcnl(nloutidx,1,nlidx,attidx)=...
             xcov(tr(tgoodidx{1}),tmod_psth(tgoodidx{1}),0,'coeff');
      end
      
      % new crap to do hinge fit
      if respfmtcode==1,
         
         % don't sum over time dimension... this returns a
         % time X 1 X latency matrix
         tmod_psth=kernpredict(hf(:,latidx,nlidx,attidx),stim,1,0,0);
         tmod_psth=squeeze(tmod_psth);
      else
         % don't rectify at all for hinge fit, regardless of nl
         % used to get hf
         tmod_psth=kernpredict(hf(:,:,nlidx,attidx),tstim,1,0);
      %elseif nlidx==1,    % no rectification
      %   tmod_psth=kernpredict(hf(:,:,nlidx,attidx),tstim,1,0);
      %elseif nlidx==2,    % rectify summed output
      %   tmod_psth=kernpredict(hf(:,:,nlidx,attidx),tstim,1,0);
      %elseif nlidx==3,    % rectify each spatial channel individually
      %   tmod_psth=kernpredict(hf(:,:,nlidx,attidx),tstim,1,0);
      end
      
      if sum(abs(tmod_psth(tgoodidx{latidx})))>0,
         hingeparms(:,latidx,nlidx,attidx)=...
             fithinge(tmod_psth(tgoodidx{latidx}),tr(tgoodidx{latidx}),0)';
         drawnow;
      else
         disp('skipping zero pred');
      end
   end
end % for attidx

xcalt=zeros(sfscount,sigcount,latcount,nlcount);
for latidx=1:latcount,
   if attcount>1,
      txc=zeros(sfscount,sigcount);
      for attidx=2:attcount,
         txc=txc+(xc(:,:,latidx,1,attidx)*tgoodlen(latidx,attidx)).^2;
      end
      xcalt(:,:,latidx,1)=sqrt(txc./sum(tgoodlen(latidx,2:end).^2));
   else
      xcalt(:,:,latidx,1)=xc(:,:,latidx,1,1);
   end
end

return


% old crap for other fitting method

figure(1);
clf
maxc=max(max(max(max(xcalt(:,:,:,:)))));
minc=min(min(min(min(xcalt(:,:,:,:)))));
sfsfit=zeros(latcount,nlcount);
sigfit=zeros(latcount,nlcount);
hf=zeros(spacecount,size(mH,2)+maxlag(1),attcount);
for latidx=1:latcount,
   for nlidx=1:nlcount,
      txc=xcalt(:,:,latidx,nlidx);
      
      %smooth txc?
      if 0,
         [XX,YY]=meshgrid(-2:2,-2:2);
         sx=1.5; sy=2;
         g=exp(-XX.^2./(2*sx.^2) - YY.^2./(2*sy.^2));
         g=g./sum(g(:));
         txc=rconv2(txc,g);
      end
      maxidx=min(find(txc==max(txc(:))));
      [sfsfit(latidx,nlidx),sigfit(latidx,nlidx)]=...
          ind2sub([sfscount,sigcount],maxidx);
      
      subplot(nlcount,latcount,(nlidx-1)*latcount+latidx);
      imagesc(xcalt(:,:,latidx,nlidx),[minc maxc]);
      
      hold on
      plot(sigfit(latidx,nlidx),sfsfit(latidx,nlidx),'x');
      hold off
      title(sprintf('FIT latidx=%d nl=%d',latidx,nlidx));
      
   end
end      
colorbar
colormap(hot);


figure(2)
showkern(hf(:,:,:,1),kernfmt,iconside);
%showkern(squeeze(mSR),kernfmt,iconside)
%showkern(squeeze(mH(:,:,end,:,:)),kernfmt,iconside)


