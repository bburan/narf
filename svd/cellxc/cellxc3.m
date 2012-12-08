% cellxc3.m
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
% 1. load and preprocess stimulus and response files as specified
%    in stimfiles,respfiles,times,stimloadcmd,stimloadparms,etc.
%    (xcloadfiles.m)
% 2. calculate stimulus-response XC and normalize. (xccore.m)
%
%
% created SVD 1/23/03 - ripped off of cellxc2.m

VCELLXC=3;

% basic stim parms. should be set early
if ~exist('stimfiles'),
   disp('no stimfile selected. can''t run cellxc.m');
   keyboard
end

% parms for XC / decorr
if ~exist('boundary'),
   boundary='zero';
end

attcount=size(times(1).start,2);  % ignore attcount>1 for the time being
attidx=1;

xcloadfiles;

if batchdata.id~=65,
   % simply do the XC. standard
   xccore;
   
else
   % experimental ARD-ish method.  iteratively remove spatial
   % dimensions that are really high in noise and re-run xcore
   
   stim0=stim;
   spacecount0=size(stim,2);
   goodidx=1:spacecount0;
   
   MAXIT=3;
   itidx=0;
   keepgoing=1;
   while keepgoing & itidx<MAXIT,
      itidx=itidx+1;
      fprintf('STARTING CORE ITERATION %d:\n',itidx);
      
      xccore;
      
      if itidx==1,
         mS0=mS;
         mSall0=mSall;
      end
      
      r=1;
      sfsidx=round(sfscount/2);
      
      smm=mH(:,:,sfsidx,r);
      sms=eH(:,:,sfsidx,r);
      smd=abs(smm)./(sms+(sms==0));  % ratio of mean to stderr
      
      smdac=mean(smd(:,2:-maxlag(1)),2);
      smdc=mean(smd(:,-maxlag(1)+1+(2:min([-maxlag(1) maxlag(2)-1]))),2);
      
      if 0,
      
      blen=size(stim,1);
      bsmean=mean(stim,1);
      tbincount=min([maxlag(2)+1 14]);
      latcount=1;
      respuse=1;
      nlcount=1;
      sfsrange=5:10:sfscount;
      sigfracs=[0.5 0.8 0.9 0.95 0.975 1.0];
      sigcount=length(sigfracs);
      mod_psth=zeros(blen,respuse,nlcount,length(sfsrange),sigcount);
      
      for resampidx=1:resampcount,
         fprintf('resampidx=%d',resampidx);
         userange=rstartidx(resampidx):rendidx(resampidx);
         if meansub,
            % subtract exploratory mean from stimulus to get DC level
            % set correctly
            tstim=stim(userange,:)-repmat(bsmean,[length(userange) 1]);
         else
            tstim=stim(userange,:);
         end
         for sfsidx=1:length(sfsrange),
            fprintf('.');
            
            for sigidx=1:sigcount,
               
               tsmdc=sort(smdc);
               goodchan=find(smdc>=tsmdc(floor(spacecount*(1-sigfracs(sigidx)))+1));
               
               thf=H(goodchan,(1:tbincount)-maxlag(1),sfsrange(sfsidx),1,resampidx);
               
               % linear prediction, rectified?
               tmod_psth=kernpredict(thf,tstim(:,goodchan)',length(goodchan),0,1);
               valididx=find(~isnan(tmod_psth(:,1)));
               tmod_psth=tmod_psth(valididx,:);
               linpred=sum(tmod_psth,2);
               
               mod_psth(userange(valididx),:,1,sfsidx,sigidx)=linpred;
               if nlcount>1,
                  thf=mean(H(goodchan,(1:tbincount)-maxlag(1),sfsrange(sfsidx),1,:),5);
                  
                  tmod_psth=kernpredict(thf,tstim(:,goodchan)',length(goodchan),0,1);
                  tmod_psth=tmod_psth(valididx,:);
                  linpred=sum(tmod_psth,2);
                  
                  mod_psth(userange(valididx),:,2,sfsidx,sigidx)=linpred;
               end
            end
         end
         fprintf('\n');
         if exist('BATQUEUEID') & BATQUEUEID>0,
            dbsetqueue(BATQUEUEID,2000+resampidx*10);
         end
      end
      
      xc=zeros(length(sfsrange),sigcount,latcount,nlcount,attcount);
      xcnl=zeros(nloutparm,latcount,nlcount,attcount);
      sfsfit=zeros(latcount,nlcount,attcount);
      sigfit=zeros(latcount,nlcount,attcount);
      for respidx=1:respuse,
         tgoodidx=find(~isnan(resp(:,respidx)));
         for nlidx=1:nlcount,
            for sfsidx=1:length(sfsrange),
               for sigidx=1:sigcount,
                  
                  xc(sfsidx,sigidx,respidx,nlidx,attidx)=...
                      xcov(mod_psth(tgoodidx,respidx,nlidx,sfsidx,sigidx),...
                           resp(tgoodidx,respidx),0,'coeff');
               end
            end
            
            % record index of maximal pred
            maxidx=min(find(xc(:,:,respidx,nlidx,attidx)==...
                            max(max(xc(:,:,respidx,nlidx,attidx)))));
            [sfsfituse,sigfituse]=ind2sub([length(sfsrange),sigcount],maxidx);
            sfsfit(respidx,nlidx,attidx)=sfsfituse;
            sigfit(respidx,nlidx,attidx)=sigfituse;
         end
      end
      
      fprintf('max xc at sigfrac=%.2f sfs=%d\n',...
              sigfracs(sigfituse(1)),sfsrange(sfsfituse(1)));
      
      goodchan=find(smdc>=tsmdc(floor(spacecount*(1-sigfracs(sigfituse(1))))+1));
      
      else
         
      % recursive to preserve original ids
      % of good spatial channels
      goodchan=find(smdc>smdac);
      
      MINRAT=0.3;
      if length(goodchan)./spacecount0 <= MINRAT,
         tsmdc=sort(smdc);
         goodchan=find(smdc>=tsmdc(ceil(spacecount-spacecount0*MINRAT)));
      end
      
      end
      
      if length(goodchan)<spacecount & itidx<MAXIT,
         goodidx=goodidx(goodchan);
         stim=stim(:,goodchan);
         spacecount=size(stim,2);
         fprintf('space reduced to %d channels.\n',spacecount);
         firstseg=1;
      else
         % didn't reduce further. stop now and proceed with other
         % fitting stuff.
         keepgoing=0;
      end
   end   
   
   
   % after final pruning iteration return STRF spatial structure to
   % orignal dimensions
   
   mS=mS0;
   mSall=mSall0;
   stim=stim0;
   clear mS0 mSall0 stim0
   
   smH=size(mH);
   mH0=mH;
   mH=zeros([spacecount0 smH(2:end)]);
   mH(goodidx,:,:,:)=mH0;
   clear mH0;
   seH=size(eH);
   eH0=eH;
   eH=zeros([spacecount0 seH(2:end)]);
   eH(goodidx,:,:,:)=eH0;
   clear eH0;
   smSR=size(mSR);
   mSR0=mSR;
   mSR=zeros([spacecount0 smSR(2:end)]);
   mSR(goodidx,:,:,:)=mSR0;
   clear mSR0
   seSR=size(eSR);
   eSR0=eSR;
   eSR=zeros([spacecount0 seSR(2:end)]);
   eSR(goodidx,:,:,:)=eSR0;
   clear eSR0
   
end

if DOCELLFIT2,
   cellfit2;
end

clear sH H ttSA tsSA2 tSR
clear tSA sSA1 sSA2 sSA0 SR

disp('cellxc3 done');
return

