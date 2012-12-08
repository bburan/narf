% xcitercore.m
%
% do repeated exclusions
%
% experimental ARD-ish method.  iteratively remove spatial
% dimensions that are really high in noise and re-run xcore
% 
% basically loops, calling xccore.m on repeatedly smaller subsets
% of spatial channels after excluding ones with low SNR
%

MAXIT=2;
itidx=0;
keepgoing=1;

%
% old .. replaced svd 5/5/04
%
goodidx=1:spacecount;

while keepgoing & itidx<MAXIT,
   
   itidx=itidx+1;
   fprintf('STARTING CORE ITERATION %d:\n',itidx);
   
   xccore;
   
   if itidx==1,
      mS0=mS;
      mSall0=mSall;
      spacecount0=spacecount;
      stim0=stim;
   end
   
   sfsidx=round(params.sfscount/2);
   
   smm=mH(:,:,sfsidx,:);
   sms=eH(:,:,sfsidx,:);
   smd=mean(abs(smm)./(sms+(sms==0)),4);  % ratio of mean to stderr

   if params.maxlag(1)<0,
      anticausalidx=2:-params.maxlag(1);
      causalidx=-params.maxlag(1)+1+...
                (2:min([-params.maxlag(1) params.maxlag(2)-1]));
      
      smdac=max(smd(:,anticausalidx),[],2);
      smdc=max(smd(:,causalidx),[],2);
      
      %smdac=abs(mean(smd(:,anticausalidx),2))+max(smd(:,anticausalidx),[],2)./2;
      %smdc=abs(mean(smd(:,causalidx),2))+max(smd(:,causalidx),[],2)./2;
      
      %smdac=abs(mean(smd(:,anticausalidx),2))+std(smd(:,anticausalidx),0,2);
      %smdc=abs(mean(smd(:,causalidx),2))+std(smd(:,causalidx),0,2);
   else
      smdc=smd;
      smdac=ones(size(smdc))*0.67;
   end
   
   % recursive to preserve original ids
   % of good spatial channels
   goodchan=find(smdc>smdac);
   
   MINRAT=0.3;
   if length(goodchan)./spacecount0 <= MINRAT,
      tsmdc=sort(smdc);
      goodchan=find(smdc>=tsmdc(ceil(spacecount-spacecount0*MINRAT)));
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

sH=size(H);
H0=H;
H=zeros([spacecount0 sH(2:end)]);
H(goodidx,:,:,:,:)=H0;
clear H0;
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

spacecount=spacecount0;

return

%
% new .. inserted svd 5/5/04
%

itidx=0;
keepgoing=1;
respiterbak=resp;

while keepgoing & itidx<MAXIT,
   itidx=itidx+1;
   fprintf('STARTING CORE ITERATION %d:\n',itidx);
   
   firstseg=1;
   xccore;
   
   r=1;
   sfsidx=round(params.sfscount/2);
   
   smm=mH(:,:,sfsidx,r);
   sms=eH(:,:,sfsidx,r);
   smd=abs(smm)./(sms+(sms==0));  % ratio of mean to stderr

   tH=shrinkage(smm(:,-params.maxlag(1)+1:end),...
                sms(:,-params.maxlag(1)+1:end),1);
   pred=kernpredict(tH,stim',1,0);
   gidx=find(~isnan(pred) & ~isnan(resp(:,1)));
   pp=pred(gidx);
   rr=resp(gidx,1);
   
   hingeparms=fithinge(pp,rr,1);
   pp=hinge(hingeparms,pp);
   
   badidx=gidx(find(abs(pp-rr)>std(pp-rr)*3));
   
   if itidx<MAXIT,
      resp(badidx,:,:,:)=nan;
      fprintf('%d outlier samples removed.\n',length(badidx));
   end
end


% restore original response

resp=respiterbak;
clear respiterbak

return



%
% old .. replaced svd 5/5/04
%
while 0 & keepgoing & itidx<MAXIT,
   
   itidx=itidx+1;
   fprintf('STARTING CORE ITERATION %d:\n',itidx);
   
   xccore;
   
   if itidx==1,
      mS0=mS;
      mSall0=mSall;
   end
   
   r=1;
   sfsidx=round(params.sfscount/2);
   
   smm=mH(:,:,sfsidx,r);
   sms=eH(:,:,sfsidx,r);
   smd=abs(smm)./(sms+(sms==0));  % ratio of mean to stderr

   if params.maxlag(1)<0,
      anticausalidx=2:-params.maxlag(1);
      causalidx=-params.maxlag(1)+1+(2:min([-params.maxlag(1) params.maxlag(2)-1]));
      
      smdac=max(smd(:,anticausalidx),[],2);
      smdc=max(smd(:,causalidx),[],2);
      
      %smdac=abs(mean(smd(:,anticausalidx),2))+max(smd(:,anticausalidx),[],2)./2;
      %smdc=abs(mean(smd(:,causalidx),2))+max(smd(:,causalidx),[],2)./2;
      
      %smdac=abs(mean(smd(:,anticausalidx),2))+std(smd(:,anticausalidx),0,2);
      %smdc=abs(mean(smd(:,causalidx),2))+std(smd(:,causalidx),0,2);
   else
      smdc=smd;
      smdac=ones(size(smdc))*0.5;
   end
   
   % recursive to preserve original ids
   % of good spatial channels
   goodchan=find(smdc>smdac);
   
   MINRAT=0.3;
   if length(goodchan)./spacecount0 <= MINRAT,
      tsmdc=sort(smdc);
      goodchan=find(smdc>=tsmdc(ceil(spacecount-spacecount0*MINRAT)));
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

spacecount=spacecount0;

return



