% xcfit2.m
%
% revised from xcfit to experiment with alternative means of selecting svd cutoff
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
% strf( ).h    :  linear kernel
% strf( ).mS   :  mean stimulus used in fit for subtracting before pred
% strf( ).name :  info about what this is
% strf( ).parms.sfsidx   : best fit idx
% strf( ).parms.sigidx   : best fit idx
% strf( ).parms.kernfmt  : kernfmt for showkern
% strf( ).parms.iconside : stim size info for showkern
% strf( ).nlparms : parms for output nl
% strf( ).nltype  : string desc of nl
% strf( ).xc      : results of fit xc test.
%
% strf array can be many dimensional (eg, nl X att)
%

disp('xcfit2.m:');

params.matchattsfs=getparm(params,'matchattsfs',0);
params.minsfsidx=getparm(params,'minsfsidx',1);
params.tbinms=getparm(params,'tbinms',14);
params.architecture=getparm(params,'architecture','linear');

fparams=params;
fparams.nrandxcov=0;

%
% figure out dimensions of relevant matrices
%

spacecount=size(H,1);
if fparams.keepneg && diff(fparams.maxlag)>0,
   tbincount=size(H,2);
elseif diff(fparams.maxlag)>0,
   tbincount=fparams.maxlag(2)+1; % min([fparams.maxlag(2)+1 14]);
else
   tbincount=size(H,2);
end
sfscount=size(H,3);
attcount=size(H,4);
resampcount=size(H,5);

if resampcount<=1,
   disp('resampcount must be >1');
   return
end

% causal time lags in the kernel
if fparams.keepneg,
   lagrange=(fparams.maxlag(1):fparams.maxlag(2))-fparams.maxlag(1)+1;
   zerobin=-fparams.maxlag(1)+1;
else
   lagrange=(1:tbincount)-fparams.maxlag(1);
   zerobin=1;
end

sfsidxtouse=(1:fparams.sfscount)';  % test all sfs cutoffs

% define a range of sigmas for combining resampled estimates (via
% shrinkage filter)
if length(fparams.sffiltsigma)>1,
   sigrange=fparams.sffiltsigma(:)';
   
elseif fparams.respfmtcode==0,
   
   %
   % USED FOR V1 ANALYSIS
   % (Adjusted down by a factor of sqrt(resampcount) after eH
   % computation was corrected to be actual standard error.
   %
   sigrange=[0 exp(linspace(log(0.2),log(1.8),fparams.sffiltsigma-1))];
else
   
   % USED FOR PFTH-DRIVEN ANALYSIS
   sigrange=exp(linspace(log(0.8),log(1.5),fparams.sffiltsigma));
end

if 0, % run fast  for debugging
   sigrange=exp([0 1]);
   sfsidxtouse=0:2:10;
   sfscount=length(sfsidxtouse);
end
sigcount=length(sigrange);

nlcount=fparams.nloutparm;
nlstr={'none','post-sum thresh','raw_nl','adapt_post','exp thresh','pre-sum thresh'};

% for PFTH preds, separate score for each latency bin
if fparams.respfmtcode==0,
   latcount=1;
else
   latcount=tbincount;
end
% based on
%predxc=ones(strfcount,attcount,latcount) * nan;
xc=zeros(sfscount,sigcount,nlcount,attcount,latcount);
nmse=zeros(sfscount,sigcount,nlcount,attcount,latcount);
mutinfo=zeros(sfscount,sigcount,nlcount,attcount,latcount);
sfsfit=zeros(latcount,nlcount,attcount);
sigfit=zeros(latcount,nlcount,attcount);
expxc=zeros(latcount,nlcount,attcount);
xcnl=zeros(latcount,nlcount,attcount);

threshparm=zeros(sfscount,sigcount,latcount,nlcount,attcount);
threshfit=zeros(latcount,nlcount,attcount);

% make sure there are no lingering fits
clear strf
clear tstrf

%samplecount=6; %number of evenly distributed sample kernels to save
%hfsample=zeros(spacecount,tbincount,samplecount,attcount);

% loop through each attentional state, sfs cutoff and sigma fudge
% factor and predict the response to the fit data
for attidx=1:attcount,
   fprintf('Fitting for attentional state %d.\n',attidx);
   
   % just use stim and resp used to fit. this requires also that rstartidx and
   % rendidx have been defined appropriately by respresmpsegs.m

   mS=mSall(:,attidx);
   tstim=stim'-repmat(mS,[1 size(stim,1)]);
   
   %
   % FIRST FIND OPTIMAL REGULARIZATION PARM
   %
   if attidx==1 | ~params.matchattsfs,
      
      fprintf('Testing fit (att/%d,sfsidx/%d):',attcount,sfscount);
      valpred=zeros(size(resp(:,attidx))).*nan;
      sigidx=1;  % test the rest later
      for sfsidx=1:params.sfscount,
         fprintf('(%d,%d)',attidx,sfsidxtouse(sfsidx));
         
         % predict the entire fit set by predicting the booted-out
         % group for each bootstrapped kernel estimate
         for resampidx=1:params.resampcount,
            
            % only use causal time lags for pred
            %thf=H(:,(1:tbincount)-fparams.maxlag(1),sfsidx,attidx,resampidx)               
            % pick either causal or all bins, depending on params.keepneg
            thf=H(:,lagrange,sfsidx,attidx,resampidx);
            
            % predict the entire fit stimulus
            tstrf=[];
            tstrf.h=thf;
            tstrf.zerobin=zerobin;
            tstrf.architecture=params.architecture;
            
            % predictions should now work for linear or order2 architectures
            tmod_psth=strfpredict(tstrf,tstim');
            
            % figure out which segments of the response were used
            % to actually fit the current resampled kernel.... 
            if resampidx==1,
               fitidx=rstartidx(2,attidx):rendidx(end,attidx);
            elseif resampidx==params.resampcount,
               fitidx=rstartidx(1,attidx):rendidx(end-1,attidx);
            else
               fitidx=[rstartidx(1,attidx):rendidx(resampidx-1,attidx) ...
                       rstartidx(resampidx+1,attidx):rendidx(end,attidx)];
            end
            
            % ... and which segments were reserved for fitting
            % hyperparameters (ie, sfsstep and sfssig)
            fitidx2=rstartidx(resampidx,attidx):rendidx(resampidx,attidx);
            
            % fit dc and gain to the estimation data
            gfitidx1=fitidx(find(~isnan(resp(fitidx,attidx)) & ...
                                 ~isnan(tmod_psth(fitidx))));
            dcgparm=fitdcgaindumb(tmod_psth(gfitidx1),resp(gfitidx1,attidx));
            
            % apply dc/g adjustment to prediction and kernel
            tmod_psth=dcgain(dcgparm,tmod_psth);
            
            H(:,:,sfsidx,attidx,resampidx)=...
                H(:,:,sfsidx,attidx,resampidx).*dcgparm(2);
            
            % save predictions of the excluded segment for choosing a
            % regularization parameter
            valpred(fitidx2)=tmod_psth(fitidx2);
         end
         
         validx=find(~isnan(resp(:,attidx)) & ~isnan(valpred));
         %xc=zeros(sfscount,sigcount,nlcount,attcount,latcount);
         xc(sfsidx,sigidx,1,attidx)=xcov(resp(validx,attidx),...
                                         valpred(validx),0,'coeff');
         
         % redefine mH and eH from re-scaled h
         mH(:,:,sfsidx,attidx)=mean(H(:,:,sfsidx,attidx,:),5);
         eH(:,:,sfsidx,attidx)=std(H(:,:,sfsidx,attidx,:),1,5) .* ...
             sqrt((params.resampcount-1)/params.resampcount .* ...
                  params.resampcount);
         
         dbsetqueue;
      end
      fprintf('\n');
   end
   
   % figure out sfsidx that gives the max xc for each latidx
   % latidx can either refer to each entry in the time lag
   % dimension (respfmtcode=1) or each response dimension
   % (respfmtcode=0). the latter is the case when you're
   % choosing the same sfsidx,sigidx for multiple kernels (eg
   % sample quantity vs pred corr analysis)
      
   % figure out max xc
   if params.matchattsfs,
      txc=xc(:,1,1,1);
   else
      txc=xc(:,1,1,attidx);
   end
   
   % option to smooth txc.  not useful any more
   if 0 & size(txc,1)>=8 & size(txc,2)>=3,
      [XX,YY]=meshgrid(-2:2,-2:2);
      sx=0.5; sy=1.0;
      g=exp(-XX.^2./(2*sx.^2) - YY.^2./(2*sy.^2));
      g=g./sum(g(:));
      txc=rconv2(txc,g);
   end
   
   if params.minsfsidx>1,
      fprintf('forced: minsfsidx=%d\n',params.minsfsidx);
   end
   
   txc(1:(params.minsfsidx-1),:)=-1;
   
   maxsfsidx=min(find(txc==max(txc(:))));
   if isempty(maxsfsidx),
      maxsfsidx=min(size(txc,1),params.minsfsidx);
   end
   sfsfit(1,:,attidx)=maxsfsidx;
   
   fprintf('optimal sfsidx=%d\n',maxsfsidx);
   sfsidx=maxsfsidx;
   
   %
   % NOW FIND OPTIMAL SHRINKAGE PARM
   %
   
   fprintf('Finding optimal shrinkage:\n');
   for sigidx=2:sigcount,
      
      % predict the entire fit set by predicting the booted-out
      % group for each bootstrapped kernel estimate
      for resampidx=1:params.resampcount,
         % current kernel is shrinkage-filtered from current (sfs,r,attidx)
         % only use causal time lags for pred!
         thf=H(:,lagrange,sfsidx,attidx,resampidx);
         
         smm=mH(:,lagrange,sfsidx,attidx);
         sms=eH(:,lagrange,sfsidx,attidx);
         smd=abs(smm)./(sms+(sms==0)) ./ sigrange(sigidx);
         if fparams.shrinkage, % old--shrinkage filter
            smd=(1-smd.^(-2));
            smd=smd.*(smd>0);
            smd(find(isnan(smd)))=0;
            thf=thf.*smd;
         else
            % new -- threshold by # of std errs
            thf=thf.*(smd>1);
         end
         if sum(thf(:))==0,
            fprintf('0');
         end
         
         % predict the entire fit stimulus
         tstrf=[];
         tstrf.h=thf;
         tstrf.zerobin=zerobin;
         tstrf.architecture=params.architecture;
         
         % predictions should now work for linear or order2 architectures
         tmod_psth=strfpredict(tstrf,tstim');

         % figure out which segments of the response were used to actually fit
         % the current resampled kernel. 
         if resampidx==1,
            fitidx=rstartidx(2,attidx):rendidx(end,attidx);
         elseif resampidx==params.resampcount,
            fitidx=rstartidx(1,attidx):rendidx(end-1,attidx);
         else
            fitidx=[rstartidx(1,attidx):rendidx(resampidx-1,attidx) ...
                    rstartidx(resampidx+1,attidx):rendidx(end,attidx)];
         end
         fitidx2=rstartidx(resampidx,attidx):rendidx(resampidx,attidx);
         
         gfitidx1=fitidx(find(~isnan(resp(fitidx,attidx)) & ...
                              ~isnan(tmod_psth(fitidx))));
         gfitidx2=fitidx2(find(~isnan(resp(fitidx2,attidx)) & ...
                               ~isnan(tmod_psth(fitidx2))));
         
         % fit dc and gain to the estimation data
         dcgparm=fitdcgaindumb(tmod_psth(gfitidx1),resp(gfitidx1,attidx));
         
         % apply dc/g adjustment to prediction and kernel
         tmod_psth=dcgain(dcgparm,tmod_psth);
         
         % save predictions of the excluded segment for choosing a regularization
         % parameter
         valpred(fitidx2)=tmod_psth(fitidx2);
      end
      
      validx=find(~isnan(resp(:,attidx)) & ~isnan(valpred));
      if length(validx)==0,
         keyboard;
      end
      
      if std(valpred(validx))>0,
         xc(sfsidx,sigidx,1,attidx)=...
             xcov(resp(validx,attidx),valpred(validx),0,'coeff');
       else
          fprintf('0!\n');
      end
       
      if exist('BATQUEUEID','var') && ~isempty(BATQUEUEID),
         % record latest (sfsidx,attidx) predicted
         dbsetqueue(BATQUEUEID);
      end
   end
   
   % figure out shrinkage parm that give max xc for each latidx
   for latidx=1:latcount,
      
      % figure out max xc
      if params.matchattsfs,
         txc=xc(maxsfsidx,:,1,1,latidx);
      else
         txc=xc(maxsfsidx,:,1,attidx,latidx);
      end
      
      maxsigidx=min(find(txc==max(txc(:))));
      if isempty(maxsigidx),
         maxsigidx=1;
      end
      sigfit(latidx,:,attidx)=maxsigidx;
   end
   
   fprintf('optimal sigidx=%d\n',maxsigidx);
   
   % generate optimal kernels and save in hf
   for nlidx=1:nlcount,
      
      % fill hf with best fit kernel according to sigfit and sfsfit
      fprintf('attidx=%d nlidx=%d: (sfs,sig)=(%d,%d)\n',...
              attidx,nlidx,sfsfit(1,nlidx,attidx),sigfit(1,nlidx,attidx));
      
      %
      % use same fit values for all latencies (in respfmtcode==1)
      %
      
      if params.respfmtcode==1,
         % only one kernel still. maybe a good idea to get rid of
         % this method and simply replace it with multiple
         % kernels of one time lag
         latloopcount=1;
         sfsfituse=floor(median(sfsfit(:,nlidx,attidx)));
         sigfituse=floor(median(sigfit(:,nlidx,attidx)));
      else
         % different kernel for each response dimension
         latloopcount=latcount;
         sfsfituse=sfsfit(:,nlidx,attidx);
         sigfituse=sigfit(:,nlidx,attidx);
         %sfsfituse=repmat(sfsfit(1,nlidx,attidx),[latcount 1]);
         %sigfituse=repmat(sigfit(1,nlidx,attidx),[latcount 1]);
      end
      xc(:,:,nlidx,:)=xc(:,:,1,:);
      
      for latidx=1:latloopcount,
         smm=mH(:,lagrange,sfsidx,attidx);
         if sigrange(sigfituse)>0 & sigfituse>1,
            sms=eH(:,lagrange,sfsidx,attidx);
            smd=abs(smm)./(sms+(sms==0)) ./ sigrange(sigfituse);
            if fparams.shrinkage, % old--shrinkage filter
               smd=(1-smd.^(-2));
               smd=smd.*(smd>0);
               smd(find(isnan(smd)))=0;
               th=smm.*smd;
            else
               % new -- threshold by # of std errs
               th=smm.*(smd>1);
            end
         else
            th=smm;
         end
         
         strf(nlidx,attidx,latidx).h=th;
         strf(nlidx,attidx,latidx).parms.sfsfit=sfsfituse(latidx);
         strf(nlidx,attidx,latidx).parms.sigfit=sigfituse(latidx);
         strf(nlidx,attidx,latidx).parms.kernfmt=fparams.kernfmt;
         strf(nlidx,attidx,latidx).parms.iconside=fparams.iconside;
         strf(nlidx,attidx,latidx).parms.tbinms=fparams.tbinms;
         strf(nlidx,attidx,latidx).parms.respfmtcode=fparams.respfmtcode;
         global DEPSTIMMAX
         if ~isempty(DEPSTIMMAX),
            strf(nlidx,attidx,latidx).parms.depstimmax=DEPSTIMMAX;
         end
         strf(nlidx,attidx,latidx).mS=mS;
         strf(nlidx,attidx,latidx).zerobin=zerobin;
         strf(nlidx,attidx,latidx).architecture=params.architecture;
         
         if exist('powunbiased','var'),
            strf(nlidx,attidx,latidx).powunbiased=...
                powunbiased(:,1,sfsfituse(latidx));
         else
            strf(nlidx,attidx,latidx).powunbiased=...
                zeros(size(th,1),1);
         end
         
         %
         % fit output nonlinearities for each h (in the
         % nlidx,latidx loop)
         %
         linpred=strfpredict(strf(nlidx,attidx,latidx),tstim');
         
         tr=resp(:,attidx);
         tgoodidx=find(~isnan(tr) & ~isnan(linpred));
         
         % fit the appropriate output NL
         switch nlstr{nlidx},
          case 'none',
            % no output NL
            tmod_psth=linpred;
            t0=min(tmod_psth(tgoodidx));
          case 'post-sum thresh',
            % rectify summed output
            tmod_psth=linpred;
            if 1 | latidx==1,
               [t0,res]=findthresh(tmod_psth(tgoodidx),tr(tgoodidx));
            else
               t0=strf(nlidx,attidx,1).nlparms;
            end
            tmod_psth=thresh(t0,tmod_psth);
          case 'raw_nl',
           % fit raw output nonlinearity.  
           disp('finding raw_nl');
           %keyboard
           if strcmp(strf(nlidx,attidx,latidx).architecture,'order2'),
              strf(nlidx,attidx,latidx)=...
                 strfdimremap(strf(nlidx,attidx,latidx),tstim',tr+mR);
              [tmod_psth,pl]=strfpredict(strf(nlidx,attidx,latidx),tstim');
              t0=strf(nlidx,attidx,latidx).nlparms;
              
              % from before strfdimremap.m :
              %strf(nlidx,attidx,latidx).architecture='stc';
              %linpred=strfpredict(strf(nlidx,attidx,latidx),tstim');
              %t0=find_raw_nl(linpred(tgoodidx,:),tr(tgoodidx),0);
              %tmod_psth=raw_nl(t0,linpred);
           else
              t0=find_raw_nl(linpred(tgoodidx,:),tr(tgoodidx)+mean(mR),0);
              tmod_psth=raw_nl(t0,linpred);
           end
           
           %xcov(double(tr(tgoodidx)),tmod_psth(tgoodidx),0,'coeff')
          case 'adapt_post',
           % temporal adaptation
           % scale linear predictions in the right ballpark
           r0=tr(tgoodidx,1)-mean(tr(tgoodidx,1));
           r1=linpred(tgoodidx)-mean(linpred(tgoodidx));
           d1=sum(r1.^2);
           if d1>0,
              scf=sum(r0.*r1)./d1;
           else
              scf=1;
           end
           
           % apply scaling constant to kernel.
           strf(nlidx,attidx,latidx).h=strf(nlidx,attidx,latidx).h .* scf;
           tmod_psth=strfpredict(strf(nlidx,attidx,latidx),tstim');
           
           %t0=findtadapt(tmod_psth(tgoodidx),tr(tgoodidx));
           %keyboard
           t0=fit_adapt_post(tmod_psth(tgoodidx),double(tr(tgoodidx)));
           
          case 'exp thresh',
            % threshold and find expansive NL exponent
            tmod_psth=linpred;
            if 1 | latidx==1,
               t0=fitexpthresh(tmod_psth(tgoodidx),double(tr(tgoodidx)),0);
            else
               t0=strf(nlidx,attidx,1).nlparms;
            end
            tmod_psth=expthresh(t0,tmod_psth);
          case 'pre-sum thresh',
           
            % rectify each spatial channel before sum
            strf(nlidx,attidx,latidx).architecture='linear_sep';
            seplinpred=strfpredict(strf(nlidx,attidx,latidx),tstim');
            tmod_psth=seplinpred;
            if 1 | latidx==1,
               t0=fithinge2(tmod_psth(tgoodidx,:),double(tr(tgoodidx)),0,4);
            else
               t0=strf(nlidx,attidx,1).nlparms;
            end
            tmod_psth=thresh(t0,tmod_psth);
         end
         
         if std(double(tr(tgoodidx)))>0 & std(tmod_psth(tgoodidx))>0,
            xcnl(latidx,nlidx,attidx)=xcov(double(tr(tgoodidx)),...
                                           tmod_psth(tgoodidx),0,'coeff');
         end
         expxc(latidx,nlidx,attidx)=xcnl(latidx,nlidx,attidx);
         
         if strcmp(nlstr{nlidx},'adapt_post')||...
               strcmp(nlstr{nlidx},'exp thresh'),
            r0=tr(tgoodidx,1)-mean(tr(tgoodidx,1));
            r1=tmod_psth(tgoodidx)-mean(tmod_psth(tgoodidx));
            d1=sum(r1.^2);
            if d1>0,
               scf=sum(r0.*r1)./d1;
            else
               scf=1;
            end
            
            % apply scaling constant to kernel.
            t0=t0 .* scf;
            strf(nlidx,attidx,latidx).h=strf(nlidx,attidx,latidx).h .* scf;
         end
         
         % same NL fit for each partial kernel
         
         strf(nlidx,attidx,latidx).nlparms=t0;
         strf(nlidx,attidx,latidx).nltype=nlstr{nlidx};
         strf(nlidx,attidx,latidx).name=...
             sprintf('%s NL: %s',basename(fparams.outfile),nlstr{nlidx});
         
      end   % for each latidx      
      
   end  % for each nlidx
end

%clear tmod_psth tgoodidx tgoodlen tstimloadparms tt txc thf
clear r0 r1 smm smd sms scf sfsidx ridx nlidx maxidx fidx
clear d1 XX YY

return


% old crap for other fitting method

figure(1);
clf
maxc=max(max(max(max(xcalt(:,:,:,:)))));
minc=min(min(min(min(xcalt(:,:,:,:)))));
sfsfit=zeros(latcount,nlcount);
sigfit=zeros(latcount,nlcount);
hf=zeros(spacecount,size(mH,2)+fparams.maxlag(1),attcount);
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


