% xcfit.m
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

disp('xcfit.m:');

params.matchattsfs=getparm(params,'matchattsfs',0);
params.minsfsidx=getparm(params,'minsfsidx',1);

fparams=params;
fparams.nrandxcov=0;

%
% figure out dimensions of relevant matrices
%

spacecount=size(mH,1);
if diff(fparams.maxlag)>0,
   tbincount=min([fparams.maxlag(2)+1 14]);
else
   tbincount=size(mH,2);
end
sfscount=size(mH,3);
respcount=size(mH,4);
attcount=size(mH,5);
lagrange=(1-fparams.maxlag(1)):(tbincount-fparams.maxlag(2)+10);

sfsidxtouse=(1:fparams.sfscount)';  % test all sfs cutoffs

% define a range of sigmas for combining resampled estimates (via
% shrinkage filter)
if length(fparams.sffiltsigma)>1,
   sigrange=fparams.sffiltsigma(:)';
   fprintf('user-specified sig range %.2f to %.2f\n',...
           sigrange(1),sigrange(end));
   
elseif fparams.respfmtcode==0,
   
   %
   % USED FOR V1 ANALYSIS
   % (Adjusted down by a factor of sqrt(resampcount) after eH
   % computation was corrected to be actual standard error.
   %
   if strcmp(params.cellid(1:3),'mod'),
      % model cell, use a wide range
      sigrange=exp(linspace(log(1.0),log(1.8),fparams.sffiltsigma));
   else
      sigrange=exp(linspace(log(0.9),log(1.8),fparams.sffiltsigma));
   end
   fprintf('using sig range %.2f to %.2f (count %d)\n',...
           sigrange(1),sigrange(end),length(sigrange));
else
   
   % USED FOR PFTH-DRIVEN ANALYSIS
   sigrange=exp(linspace(log(0.9),log(1.8),fparams.sffiltsigma));
end

if 0, % run fast  for debugging
   sigrange=exp([0 1]);
   sfsidxtouse=0:2:10;
   sfscount=length(sfsidxtouse);
end
sigcount=length(sigrange);

nlcount=fparams.nloutparm;
nlstr={'none','post-sum thresh','adapt_post','pre-sum thresh','exp thresh','sigmoid'};

% for PFTH preds, separate score for each latency bin
if fparams.respfmtcode==0,
   latcount=respcount;
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

clear tstrf

%samplecount=6; %number of evenly distributed sample kernels to save
%hfsample=zeros(spacecount,tbincount,samplecount,attcount);

% loop through each attentional state, sfs cutoff and sigma fudge
% factor and predict the response to the fit data
if params.respfmtcode==0,
   MAXPREDLEN=10000;
elseif params.respfmtcode==1,
   MAXPREDLEN=100000;
end

for attidx=1:attcount,
   fprintf('Fitting for attentional state %d.\n',attidx);
   
   
   if exist('fdata','var'),
      
      % do nothing
      
   elseif ~exist('VCELLXC','var') | VCELLXC==1 | fparams.fitfrac>0,
      % load stim and response segments to fit. reload the stimulus
      % and response for each attentional state
      
      [fdata.stim,fdata.resp]=xcloadstimresp(fitfile(:,attidx),...
                                             fitstartframe(:,attidx),...
                                             fitstopframe(:,attidx),...
                                             fparams);
   elseif size(resp,1)>MAXPREDLEN,
      rgoodidx=find(~isnan(resp(:,1)));
      rminidx=rgoodidx(1);
      if length(rgoodidx)>=MAXPREDLEN,
         rmaxidx=rgoodidx(MAXPREDLEN);
      else
         rmaxidx=rgoodidx(end);
      end
      fprintf('trimming exp data to frames %d-%d for fitting\n',...
              rminidx,rmaxidx);
      fdata.resp=resp(rminidx:rmaxidx,:);
      fdata.stim=stim(rminidx:rmaxidx,:);
   else
      fdata.stim=stim;
      fdata.resp=resp;
   end
   
   fprintf('Testing fit (att/%d,sfsidx/%d):',attcount,sfscount);
   respsave=fdata.resp;
   for sfsidx=params.minsfsidx:sfscount,
      fprintf('(%d,%d)',attidx,sfsidxtouse(sfsidx));
      
      for respidx=1:respcount,
         for sigidx=1:sigcount,
            % current kernel is shrinkage-filtered from current (sfs,r,attidx)
            % only use causal time lags for pred!
            if fparams.resampcount>1,
               valtimes=(1:tbincount)-fparams.maxlag(1);
               thf=shrinkage(mH(:,valtimes,sfsidx,respidx,attidx),...
                             eH(:,valtimes,sfsidx,respidx,attidx),...
                             sigrange(sigidx),fparams.shrinkage);
               
               %smm=mH(:,(1:tbincount)-fparams.maxlag(1),sfsidx,respidx,attidx);
               %sms=eH(:,(1:tbincount)-fparams.maxlag(1),sfsidx,respidx,attidx);
               %smd=abs(smm)./(sms+(sms==0)) ./ sigrange(sigidx);
               %if fparams.shrinkage, % old--shrinkage filter
               %   smd=(1-smd.^(-2));
               %   smd=smd.*(smd>0);
               %   smd(find(isnan(smd)))=0;
               %   thf=smm.*smd;
               %else
               %   % new -- threshold by # of std errs
               %   thf=smm.*(smd>1);
               %end
            else
               % sigidx has no influence. assuming sigcount=1
               thf=mH(:,(1:tbincount)-fparams.maxlag(1),sfsidx,respidx,attidx);
            end
            
            for nlidx=1:nlcount,
               tstrf(nlidx,sigidx).h=thf;
               tstrf(nlidx,sigidx).nlparms=[0; 0];
               tstrf(nlidx,sigidx).nltype=nlstr{nlidx};
               tstrf(nlidx,sigidx).mS=mSall(:,respidx,attidx);
               tstrf(nlidx,sigidx).parms.kernfmt=fparams.kernfmt;
               tstrf(nlidx,sigidx).parms.iconside=fparams.iconside;
               tstrf(nlidx,sigidx).name=...
                   sprintf('%s sfs=%d sig=%d nl=%d',fparams.cellid,...
                           sfsidx,sigidx,nlidx);
            end
         end
         
         % temporarily turned off
         if params.matchattsfs,
            fdata.resp=respsave(:,1);
         else
            fdata.resp=respsave(:,attidx);
         end
         
         tfitres=xcval(tstrf(:),fparams,fdata);
         
         xc(sfsidx,:,:,attidx,respidx)= ...
             permute(reshape(tfitres.predxc,1,nlcount,sigcount,1),[1 3 2 4]);
         nmse(sfsidx,:,:,attidx,respidx)= ...
             permute(reshape(tfitres.predmse,1,nlcount,sigcount,1),[1 3 2 4]);
      end
      
   end     % for sfsidx
   
   if not(isempty(BATQUEUEID)),
      % record latest (sfsidx,attidx) predicted
      dbsetqueue(BATQUEUEID,1);
   end
   
   fprintf('\n');
   fdata.resp=respsave;
   
   %xc=zeros(sfscount,sigcount,nlcount,attcount,latcount);
   
   % generate optimal kernels and save in hf
   for nlidx=1:nlcount,
      
      % figure out max xc for each latidx
      % latidx can either refer to each entry in the time lag
      % dimension (respfmtcode=1) or each response dimension
      % (respfmtcode=0). the latter is the case when you're
      % choosing the same sfsidx,sigidx for multiple kernels (eg
      % sample quantity vs pred corr analysis)
      
      for latidx=1:latcount,
         
         % figure out max xc
         if 0 & params.matchattsfs,
            txc=xc(:,:,nlidx,1,latidx);
         elseif 0,
            txc=xc(:,:,nlidx,1,latidx);
         else
            txc=xc(:,:,nlidx,attidx,latidx);
         end
         
         if 0,
            mtxc=mean(txc(:,2:end-1),2);
            maxsfsidx=min(find(mtxc==max(mtxc)));
            maxsigidx=max(find(txc(maxsfsidx,:)==max(txc(maxsfsidx,:))));
         else
            
            % smooth txc
            if size(txc,1)>=8 & size(txc,2)>=3,
               [XX,YY]=meshgrid(-1:1,-2:2);
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
         
         
         if nlidx>=5,
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
      
      if params.respfmtcode==1
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
      
      for latidx=1:latloopcount,
         
         if params.respfmtcode==1,
            if fparams.resampcount>1,
               smm=mH(:,:,sfsfituse(latidx),1,attidx);
               sms=eH(:,:,sfsfituse(latidx),1,attidx) .* ...
                   sigrange(sigfituse(latidx));
               smd=abs(smm)./(sms+(sms==0));
               if fparams.shrinkage,
                  smd=1-smd.^(-2);
                  smd=smd.*(smd>0);
                  smd(find(isnan(smd)))=0;
                  th=smm.*smd;
               else
                  th=smm.*(smd>1);
               end
            else
               th=mH(:,:,sfsfituse(latidx),1,attidx);
            end
            
         else
            if fparams.resampcount>1,
               valtimes=(1:tbincount)-fparams.maxlag(1);
               th=shrinkage(mH(:,valtimes,sfsfituse(latidx),latidx,attidx),...
                            eH(:,valtimes,sfsfituse(latidx),latidx,attidx),...
                            sigrange(sigfituse(latidx)));
            else
               th=mH(:,(1:tbincount)-fparams.maxlag(1),...
                    sfsfituse(latidx),latidx,attidx);
            end
         end
         
         strf(nlidx,attidx,latidx).h=th;
         strf(nlidx,attidx,latidx).parms.sfsfit=sfsfituse(latidx);
         if exist('pctvar','var'),
            strf(nlidx,attidx,latidx).parms.pctvar=pctvar(sfsfituse(latidx));
         else
            strf(nlidx,attidx,latidx).parms.pctvar=0;
         end
         strf(nlidx,attidx,latidx).parms.sigfit=sigfituse(latidx);
         strf(nlidx,attidx,latidx).parms.kernfmt=fparams.kernfmt;
         strf(nlidx,attidx,latidx).parms.iconside=fparams.iconside;
         strf(nlidx,attidx,latidx).parms.tbinms=fparams.tbinms;
         strf(nlidx,attidx,latidx).parms.respfmtcode=fparams.respfmtcode;
         strf(nlidx,attidx,latidx).mS=mSall(:,latidx,attidx);
         if exist('powunbiased','var'),
            strf(nlidx,attidx,latidx).powunbiased=...
                powunbiased(:,1,sfsfituse(latidx));
         else
            strf(nlidx,attidx,latidx).powunbiased=...
                zeros(size(th,1));
         end
         
         %
         % fit output nonlinearities for each h (in the
         % nlidx,latidx loop)
         %
         
         % subtract appropriate mean
         if fparams.meansub,
            tstim=fdata.stim'-repmat(strf(nlidx,attidx,latidx).mS,...
                                     [1 size(fdata.stim,1)]);
         else
            tstim=fdata.stim';
         end
         tr=fdata.resp(:,attidx);
         
         % scale linear kernel to give predictions in the right ballpark
         seplinpred=kernpredict(strf(nlidx,attidx,latidx).h,...
                                tstim,spacecount,0);
         linpred=sum(seplinpred,2);
         tgoodidx=find(~isnan(tr) & ~isnan(linpred));
         
         if nlidx==3,
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
            linpred=kernpredict(strf(nlidx,attidx,latidx).h,tstim,1,0);
         end
         
         % fit the appropriate output NL
         if nlidx==1,
            % no output NL
            tmod_psth=linpred;
            t0=min(tmod_psth(tgoodidx));
         elseif nlidx==2,
            % rectify summed output
            tmod_psth=linpred;
            if 1 | latidx==1,
               [t0,res]=findthresh(tmod_psth(tgoodidx),tr(tgoodidx));
            else
               t0=strf(nlidx,attidx,1).nlparms;
            end
            tmod_psth=thresh(t0,tmod_psth);
         elseif nlidx==3,
            % sigmoid
            tmod_psth=linpred;
            t0=fit_adapt_post(tmod_psth(tgoodidx),tr(tgoodidx));
            tmod_psth=adapt_post(t0,tmod_psth);
         elseif nlidx==4,
            % rectify each spatial channel
            tmod_psth=seplinpred;
            if 1 | latidx==1,
               t0=fithinge2(tmod_psth(tgoodidx,:),tr(tgoodidx),0,4);
            else
               t0=strf(nlidx,attidx,1).nlparms;
            end
            tmod_psth=thresh(t0,tmod_psth);
         elseif nlidx==5,
            % threshold and find expansive NL exponent
            tmod_psth=linpred;
            if 1 | latidx==1,
               t0=fitexpthresh(tmod_psth(tgoodidx),tr(tgoodidx),0);
            else
               t0=strf(nlidx,attidx,1).nlparms;
            end
            tmod_psth=expthresh(t0,tmod_psth);
         elseif nlidx==6,
            % sigmoid
            tmod_psth=linpred;
            if 1 | latidx==1,
               t0=fitsigmoid(tmod_psth(tgoodidx),tr(tgoodidx),0);
            else
               t0=strf(nlidx,attidx,1).nlparms;
            end
            tmod_psth=sigmoid(t0,tmod_psth);
         end
         
         if std(tr(tgoodidx))>0 & std(tmod_psth(tgoodidx))>0,
            xcnl(latidx,nlidx,attidx)=xcov(tr(tgoodidx),...
                                           tmod_psth(tgoodidx),0,'coeff');
         end
         expxc(latidx,nlidx,attidx)=xcnl(latidx,nlidx,attidx);
         
         if nlidx~=3 & nlidx~=4 & nlidx~=6,
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
end % for attidx

clear tresp tmod_psth tgoodidx tgoodlen tstimloadparms tt txc thf
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


