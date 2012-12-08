% cellfit.m
%
% we have two things to fit: 
% 1. the svd regularization parameter (sfsidx)
% 2. the shrinkage filter fudge factor (sfiltsigma)
%
% specs (so this is interchangable with other fit routines):
% inputs: mH, eH (mean and std error of kernels)
%         params (bunch of crap listed in cellxcnodb.m)
%         stim,resp (if fitfrac=0, exploratory fit data, used for
%                   fitting other parms as well. otherwise, data is
%                   loaded using params.times(2))
%
% strf structure:
% strf( ).h     :  linear kernel
% strf( ).name  :  info about what this is
% strf( ).parms.sfsidx  : best fit idx
% strf( ).parms.sigidx  : best fit idx
% strf( ).nlparms : parms for output nl
% strf( ).nltype : string desc of nl
% strf( ).xc : results of fit xc test.
%
% strf array can be many dimensional (eg, nl X att)
%

disp('cellfit.m:');

keyboard

%
% figure out dimensions of relevant matrices
%
spacecount=size(mH,1);
if diff(params.maxlag)>0,
   tbincount=min([params.maxlag(2)+1 14]);
else
   tbincount=size(mH,2);
end
sfscount=size(mH,3);
respcount=size(mH,4);
attcount=size(mH,5);
lagrange=(1-params.maxlag(1)):(tbincount-params.maxlag(2)+10);

sfsidxtouse=(1:params.sfscount)';  % test all sfs cutoffs

% define a range of sigmas for combining resampled estimates (via
% shrinkage filter)
if params.respfmtcode==0,
   if params.sffiltsigma==1,
      sigrange=4;
   else
      %
      % USED FOR V1 ANALYSIS
      % (Adjusted down by a factor of sqrt(resampcount) after eH
      % computation was corrected to be actual standard error.
      %
      sigrange=exp(linspace(log(1.0),log(2.0),params.sffiltsigma));
   end
   
else   
   % USED FOR PFTH-DRIVEN ANALYSIS
   sigrange=exp(linspace(0,log(1.5),params.sffiltsigma));
end

if 0, % run fast  for debugging
   sigrange=exp([0 1]);
   sfsidxtouse=0:2:10;
   sfscount=length(sfsidxtouse);
end
sigcount=length(sigrange);

nlcount=params.nloutparm;
nlstr={'none','post-sum thresh','exp thresh','pre-sum thresh'};

% for PFTH preds, separate score for each latency bin
if params.respfmtcode==0,
   latcount=1;
else
   latcount=tbincount;
end
nmse=zeros(sfscount,sigcount,latcount,nlcount,attcount);
mutinfo=zeros(sfscount,sigcount,latcount,nlcount,attcount);
hingeparms=zeros(3,latcount,nlcount,attcount);
tgoodlen=zeros(latcount,attcount);
sfsfit=zeros(latcount,nlcount,attcount);
sigfit=zeros(latcount,nlcount,attcount);
expxc=zeros(latcount,nlcount,attcount);
xcnl=zeros(latcount,nlcount,attcount);
%hf=zeros(spacecount,tbincount,nlcount,attcount);
%xc=zeros(sfscount,sigcount,latcount,nlcount,attcount);
threshparm=zeros(sfscount,sigcount,latcount,nlcount,attcount);
threshfit=zeros(latcount,nlcount,attcount);

%samplecount=6; %number of evenly distributed sample kernels to save
%hfsample=zeros(spacecount,tbincount,samplecount,attcount);

for attidx=1:attcount,
   fprintf('Fitting for attentional state %d.\n',attidx);
   
   if ~exist('VCELLXC') | VCELLXC==1 | params.fitfrac>0,
      % load stim and response segments to fit. reload the stimulus
      % and response for each attentional state
      
      [stim,resp]=xcloadstimresp(fitfile(:,attidx),...
                                 fitstartframe(:,attidx),...
                                 fitstopframe(:,attidx),...
                                 params);
   else
      if length(resp)>5000,
         disp('trimming fit data to 5000 frames');
         resp=resp(1:5000,:);
         stim=stim(1:5000,:);
      end
   end
   
   
   
   
   % figure out actual valid range of response being predicted
   tr=resp(:,:,attidx);
   
   % smooth actual psth if predsmoothsigma > 0
   if params.predsmoothsigma > 0 & params.respfmtcode==0,
      %& cellfiledata(fitfile).repcount==1,
      %fprintf('Smoothing actual response with optimal filter...\n');
      %pfilt=[1/9 2/9 1/3 2/9 1/9]';
      
      fprintf('Smoothing fit response sigma=%.1f bins...\n',...
              params.predsmoothsigma);
      tt=(-10:10)';
      pfilt=exp(-tt.^2/(2*params.predsmoothsigma.^2))./...
            (sqrt(2*pi)*params.predsmoothsigma);
      pfilt=pfilt(find(pfilt > 0.1*max(pfilt)));
      fr=tr;
      fr(isnan(fr))=0;
      fr=conv2(fr,pfilt,'same');
      tr(~isnan(tr))=fr(~isnan(tr));
   end
   
   % normalize valid bins of actual response
   tgoodidx={};
   for r=1:latcount,
      tgoodidx{r}=find(~isnan(tr(:,r)));
      tgoodlen(r,attidx)=length(tgoodidx{r});
   end
   
   % loop through each attentional state, sfs cutoff and sigma fudge
   % factor and predict the response to the fit data
   for r=1:respcount,
      if params.meansub,
         % subtract exploratory mean from stimulus to get DC level
         % set correctly
         tstim=stim'-repmat(mSall(:,r,attidx),[1 size(stim,1)]);
      else
         tstim=stim';
      end
      for sfsidx=1:sfscount,
         fprintf('Testing fit (att=%d r=%d sfsidx=%d sig=',...
                 attidx,r,sfsidxtouse(sfsidx));
         
         for sigidx=1:sigcount,
            fprintf('%.2f ',sigrange(sigidx));
            
            % current kernel is shrinkage-filtered from current (sfs,r,attidx)
            % only use causal time lags for pred!
            if params.resampcount>1,
               smm=mH(:,(1:tbincount)-params.maxlag(1),sfsidx,r,attidx);
               sms=eH(:,(1:tbincount)-params.maxlag(1),sfsidx,r,attidx);
               smd=abs(smm)./(sms+(sms==0)) ./ sigrange(sigidx);
               if params.shrinkage, % old--shrinkage filter
                  smd=(1-smd.^(-2));
                  smd=smd.*(smd>0);
                  smd(find(isnan(smd)))=0;
                  thf=smm.*smd;
               else
                  % new -- threshold by # of std errs
                  thf=smm.*(smd>1);
               end
            else
               % sigidx has no influence. assuming sigcount=1
               thf=mH(:,(1:tbincount)-params.maxlag(1),sfsidx,r,attidx);
            end
            
            for nlidx=1:nlcount,
               if params.respfmtcode==1,
                  tmod_psth=kernpredict(thf,tstim,1,0,0);
                  if nlidx==2,
                     tmod_psth(tmod_psth<0)=0;
                  end
                  %tmod_psth=kernpredict(thf,tstim,spacecount,1,0);
                  tmod_psth=squeeze(tmod_psth);
                  
               elseif nlidx==1,    % no rectification
                  seplinpred=kernpredict(thf,tstim,spacecount,0);
                  linpred=sum(seplinpred,2);
                  tmod_psth=linpred;
                  
               elseif nlidx==2,    % rectify summed output
                  %tmod_psth=kernpredict(thf,tstim,1,1);
                  % to speed things up, don't re-predict. just
                  % rectify previously generated linear prediction
                  tmod_psth=linpred.*(linpred>0);
                  
               elseif nlidx==3,    % no rectification
                  % linpred has already been calc'd for nlidx=1
                  tmod_psth=linpred;
                  
               elseif nlidx==4,    % rectify each spatial channel individually
                  %tmod_psth=kernpredict(thf,tstim,spacecount,1);
                  tmod_psth=seplinpred;
                  tmod_psth(find(tmod_psth<0))=0;
                  tmod_psth=sum(tmod_psth,2);
                  
               elseif 0,
                  %optionally: fit NL on each kernel!
                  tmod_psth=linpred;
                  [t0,res]=findthresh(tmod_psth(tgoodidx{1}),tr(tgoodidx{1}));
                  threshparm(sfsidx,sigidx,1,nlidx,attidx)=t0;
                  tmod_psth(tgoodidx{1})=thresh(t0,tmod_psth(tgoodidx{1}));
               end
               
               % remove leading bins that don't matter
               for ridx=1:latcount,
                  % normalize for appropriate correlation
                  % coefficient calculation
                  pp=tmod_psth(tgoodidx{ridx},ridx);
                  rr=tr(tgoodidx{ridx},ridx);
                  if std(pp)>0 & std(rr)>0,
                     xc(sfsidx,sigidx,ridx,nlidx,attidx)=...
                         xcov(rr,pp,0,'coeff');
                     rr=rr-mean(rr);
                     pp=(pp-mean(pp))./std(pp).*std(rr);
                     nmse(sfsidx,sigidx,ridx,nlidx,attidx)=...
                         var(pp-rr)./var(rr);
                     
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
         fprintf(')\n');
         
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
         
         if 0,
            
            mtxc=mean(txc(:,2:end-1),2);
            maxsfsidx=min(find(mtxc==max(mtxc)));
            maxsigidx=max(find(txc(maxsfsidx,:)==max(txc(maxsfsidx,:))));
            
         else
            
            %smooth txc?
            if size(txc,1)>=8 & size(txc,2)>=3,
               [XX,YY]=meshgrid(-2:2,-2:2);
               sx=0.5; sy=1.0;
               % fvvs
               % sx=1.5; sy=2.0;
               g=exp(-XX.^2./(2*sx.^2) - YY.^2./(2*sy.^2));
               g=g./sum(g(:));
               txc=rconv2(txc,g);
            end
            
            maxidx=min(find(txc==max(txc(:))));
            [maxsfsidx,maxsigidx]=ind2sub([sfscount,sigcount],maxidx);
         end
         
         sfsfit(latidx,nlidx,attidx)=maxsfsidx;
         sigfit(latidx,nlidx,attidx)=maxsigidx;
         strf(nlidx,attidx).parms.sfsfit(latidx)=maxsfsidx;
         strf(nlidx,attidx).parms.sigfit(latidx)=maxsigidx;
         
         % PSTH or PFTH model?
         strf(nlidx,attidx).respfmtcode=params.respfmtcode;
         
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
      
      if latcount>1,
         strf(nlidx,attidx).h=zeros(spacecount,tbincount);
         for latidx=1:latcount,
            if params.resampcount>1,
               smm=mH(:,latidx,sfsfituse,1,attidx);
               sms=eH(:,latidx,sfsfituse,1,attidx) .* sigrange(sigfituse);
               smd=abs(smm)./(sms+(sms==0));
               if params.shrinkage,
                  smd=1-smd.^(-2);
                  smd=smd.*(smd>0);
                  smd(find(isnan(smd)))=0;
                  strf(nlidx,attidx).h(:,latidx)=smm.*smd;
               else
                  strf(nlidx,attidx).h(:,latidx)=smm.*(smd>1);
               end
            else
               strf(nlidx,attidx).h(:,latidx)=mH(:,latidx,sfsfituse,1,attidx);
            end
         end
      else
         if params.resampcount>1,
            smm=mH(:,(1:tbincount)-params.maxlag(1),sfsfituse,1,attidx);
            sms=eH(:,(1:tbincount)-params.maxlag(1),sfsfituse,1,attidx);
            smd=abs(smm)./(sms+(sms==0)) ./ sigrange(sigfituse);
            if params.shrinkage,
               smd=1-smd.^(-2);
               smd=smd.*(smd>0);
               smd(find(isnan(smd)))=0;
               strf(nlidx,attidx).h=smm.*smd;
            else
               strf(nlidx,attidx).h=smm.*(smd>1);
            end
         else
            strf(nlidx,attidx).h=...
                mH(:,(1:tbincount)-params.maxlag(1),sfsfituse,1,attidx);
         end
      end
      
      % scale linear kernel to give predictions in the right ballpark
      seplinpred=kernpredict(strf(nlidx,attidx).h,tstim,spacecount,0);
      tmod_psth=sum(seplinpred,2);
      
      r0=tr(tgoodidx{1},1)-mean(tr(tgoodidx{1},1));
      r1=tmod_psth(tgoodidx{1})-mean(tmod_psth(tgoodidx{1}));
      d1=sum(r1.^2);
      if d1>0,
         scf=sum(r0.*r1)./d1;
      else
         scf=1;
      end
      
      % scale kernel to minimize absolute error
      strf(nlidx,attidx).h=scf.*strf(nlidx,attidx).h;
      strf(nlidx,attidx).mS=mSall(:,1,attidx);
      
      keyboard
      
      % fit output nonlinearities for each hf
      if nlidx==1,
         % no rectification
         seplinpred=kernpredict(strf(nlidx,attidx).h,tstim,spacecount,0);
         linpred=sum(seplinpred,2);
         tmod_psth=linpred;
         t0=min(tmod_psth(:));
      elseif nlidx==2,
         % rectify summed output
         tmod_psth=linpred;
         t0=findthresh(tmod_psth(tgoodidx{1}),tr(tgoodidx{1}),0);
      elseif nlidx==3,
         % rectify each spatial channel
         tmod_psth=linpred;
         [t0,res]=fitexpthresh(tmod_psth(tgoodidx{1},:),...
                             tr(tgoodidx{1}));
      elseif nlidx==4,
         % rectify each spatial channel
         tmod_psth=seplinpred;
         [t0,res]=findthresh(tmod_psth(tgoodidx{1},:),...
                             tr(tgoodidx{1}));
      end
      
      tmod_psth=thresh(t0,tmod_psth);
      threshfit(1:lenth(t0),nlidx,attidx)=t0;
      xcnl(latidx,nlidx,attidx)=xcov(tr(tgoodidx{1}),...
                                     tmod_psth(tgoodidx{1}),0,'coeff');
      expxc(latidx,nlidx,attidx)=xcnl(latidx,nlidx,attidx);
      
      strf(nlidx,attidx).nlparms=t0;
      strf(nlidx,attidx).nltype=nlstr{nlidx};
      strf(nlidx,attidx).name=...
          sprintf('%s NL: %s',basename(params.outfile),nlstr{nlidx});
      
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

clear tresp tmod_psth tgoodidx tgoodlen tstimloadparms tt txc thf
clear r0 r1 smm smd sms stim scf sfsidx ridx nlidx maxidx fidx
clear d1 XX YY 

return


% old crap for other fitting method

figure(1);
clf
maxc=max(max(max(max(xcalt(:,:,:,:)))));
minc=min(min(min(min(xcalt(:,:,:,:)))));
sfsfit=zeros(latcount,nlcount);
sigfit=zeros(latcount,nlcount);
hf=zeros(spacecount,size(mH,2)+params.maxlag(1),attcount);
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


