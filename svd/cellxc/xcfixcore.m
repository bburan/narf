% xcfixcore3.m  : fixation-based rc
%
% created SVD 10/18/03
%

VCELLXC=3;

% optionally smooth pfths before fitting variable gain/dc terms
TRIMFRAMES=0;
eresp=respbak(:,(TRIMFRAMES+1):end);
vresp=cdata.respbak(:,(TRIMFRAMES+1):end);
fixlen=size(eresp,2);

if params.smoothtime,
   realvalcount=nansum((eresp')~=round(eresp'));
   realvalidx=find(realvalcount>0 & ~isnan(realvalcount));
   intvalidx=find(realvalcount==0);
   
   if length(intvalidx>0),
      fprintf('smoothing %d single trials of fit data a lot\n',...
              length(intvalidx));
      g=[.05 .2 .5 .2 .05];
      %g=[.1 .2 .4 .2 .1];
      eresp(intvalidx,:)=conv2(eresp(intvalidx,:),g,'same');
   end
   if length(realvalidx>0),
      fprintf('smoothing %d multi-trials of fit data less\n',...
              length(realvalidx));
      g=[.15 .7 .15];
      %g=[.2 .6 .2];
      eresp(realvalidx,:)=conv2(eresp(realvalidx,:),g,'same');
   end
   eresp(:,1)=respbak(:,TRIMFRAMES+1);
   eresp(:,end)=respbak(:,end);
else
   disp('not smoothing in time');
end

% initally fit only resp summed over the full fixation
% (no trans/sust since that just seems to be noisey)

resp=aresp';
rsize=size(resp);
resplen=rsize(1);
respcount=rsize(2);
params.maxlag=[-2 1];
resp0idx=-params.maxlag(1)+1;
framecount=diff(params.maxlag)+1;

cdata.resp=varesp';

nnidx=find(~isnan(aresp'));
vnnidx=find(~isnan(varesp'));


% set parms to get xccore to run correctly
spacecount=size(stim,2);
firstseg=1;

params.smoothtime=0;
params.decorrtime=1;

xccore;

%
% PICK OPTIMAL CUTOFF BY BRUTE FORCE
%

tnloutparm=params.nloutparm;
params.nloutparm=1;  % only do linear fit here

xcfit2;
%xcfit;

% use regularization parameters selected by xcfit2
showidx1=strf(1).parms.sfsfit;
showsig1=strf(1).parms.sigfit;
sigtouse=sigrange(strf(1).parms.sigfit);

predpower=xc(:,:,1);
fprintf('from xcfit, maxsfsidx=%d, sigtouse=%.1f\n',showidx1,sigtouse);

params.nloutparm=tnloutparm;
predpower2=zeros(params.sfscount,params.resampcount);
vpredpower2=zeros(params.sfscount,params.resampcount);
showidx2=1;
showidx3=1;

% crunch down resampled kernels, weigh by shrinkage
thf=shrinkage(mH,eH,sigtouse);

% different output nonlinearities
nlstr={'one var gain','causal 2fix','anticausal 2fix','pn gain','p-only gain'};
nlcount=min([length(nlstr) params.nloutparm]);


%figure(2);
%showkern(thf(:,:,round(linspace(1,params.sfscount,5))),params.kernfmt)

% matrices for saving predictions and correlation scores for each sfs
rpred=zeros(resplen,framecount,params.sfscount);
posrpred=zeros(resplen,framecount,params.sfscount);
negrpred=zeros(resplen,framecount,params.sfscount);
fakeposrpred=zeros(resplen,framecount,params.sfscount);
fakenegrpred=zeros(resplen,framecount,params.sfscount);
vpred=zeros(vallen,framecount,params.sfscount,nlcount);
posvpred=zeros(vallen,framecount,params.sfscount);
negvpred=zeros(vallen,framecount,params.sfscount);
fakeposvpred=zeros(vallen,framecount,params.sfscount);
fakenegvpred=zeros(vallen,framecount,params.sfscount);
xc=zeros(1,params.sfscount,nlcount);
nlparm=zeros(1,params.sfscount,nlcount);
xcth=zeros(1,params.sfscount,nlcount);
allpredxc=zeros(1,params.sfscount,nlcount);
predxc=zeros(1,nlcount);
predp=zeros(1,nlcount);
prederr=zeros(1,nlcount);
predfix=zeros(1,nlcount);
predfixerr=zeros(1,nlcount);
sfsrec=zeros(nlcount,1);

laggain=zeros(fixlen,4,params.sfscount,nlcount);
allvpred=zeros(length(vnnidx),fixlen,params.sfscount,nlcount);
allrpred=zeros(length(nnidx),fixlen,params.sfscount,nlcount);

% predict est and val responses with each regularized kernel
fprintf('fitting nls:\n  sfsidx=');
for ii=1:params.sfscount,
   
   fprintf('%d ',ii);
   
   stim0=stim-repmat(mSall',[resplen 1]);
   rpred(:,:,ii)=stim0*thf(:,:,ii);
   
   cstim0=cdata.stim-repmat(mSall',[size(cdata.stim,1) 1]);
   vpred(:,:,ii)=cstim0*thf(:,:,ii);
   
   % shift pred responses in time to make fixation lags line up
   % appropriately
   for respidx=1:framecount,
      ridxdiff=resp0idx-respidx;
      
      if respidx<resp0idx,
         r1=[rpred((ridxdiff+1):end,respidx,ii); ones(ridxdiff,1).*nan];
      else
         r1=[ones(-ridxdiff,1).*nan; rpred(1:(end+ridxdiff),respidx,ii)];
      end
      rpred(:,respidx,ii)=r1;
      
      if respidx<resp0idx,
         r1=[vpred((ridxdiff+1):end,respidx,ii); ones(ridxdiff,1).*nan];
      else
         r1=[ones(-ridxdiff,1).*nan; vpred(1:(end+ridxdiff),respidx,ii)];
      end
      vpred(:,respidx,ii)=r1;
   end
   
   if 1,
      tthf=thf(:,:,ii);
      posrpred(:,:,ii)=stim0*(tthf.*(tthf>0));
      negrpred(:,:,ii)=stim0*(tthf.*(tthf<0));
      
      posvpred(:,:,ii)=cstim0*(tthf.*(tthf>0));
      negvpred(:,:,ii)=cstim0*(tthf.*(tthf<0));
      
      ms=2:2:spacecount;
      fakeposrpred(:,:,ii)=stim0(:,ms)*tthf(ms,:);
      fakenegrpred(:,:,ii)=stim0(:,ms-1)*tthf(ms-1,:);
      
      fakeposvpred(:,:,ii)=cstim0(:,ms)*tthf(ms,:);
      fakenegvpred(:,:,ii)=cstim0(:,ms-1)*tthf(ms-1,:);
      
   else
      ms=round(spacecount/2);
      posrpred(:,:,ii)=stim0(:,1:ms)*thf(1:ms,:,ii);
      negrpred(:,:,ii)=stim0(:,(1:ms)+ms)*thf((1:ms)+ms,:,ii);
      
      posvpred(:,:,ii)=cstim0(:,1:ms)*thf(1:ms,:,ii);
      negvpred(:,:,ii)=cstim0(:,(1:ms)+ms)*thf((1:ms)+ms,:,ii);
      
   end
   
   % fit each output kernel and adjust validation preds accordingly
   for nlidx=1:nlcount,
      
      respidx=resp0idx;
      
      if sum(abs(thf(:,respidx,ii)))>0,
         fullrpred=repmat(rpred(nnidx,respidx,ii),[1 fixlen]);
         fullvpred=repmat(vpred(vnnidx,respidx,ii),[1 fixlen]);
         
         tr=eresp(nnidx,:);
         
         switch nlidx,
            
          case 1,
           % variable gain, single threshold
           % (for linear, no short time don't do anything)
           
           for lagidx=1:fixlen,
              lagnnidx=find(~isnan(eresp(:,lagidx)));
              
              r0=eresp(lagnnidx,lagidx);
              r1=rpred(lagnnidx,respidx,ii);
              
              % no restrictions on threshold. float for each lagidx
              if sum(abs(r1))==0,
                 nlparms=[0 0];
              else
                 nlparms=fithinge3(r1,r0,0,1);
              end
              
              fullrpred(:,lagidx)=hinge3(nlparms,rpred(nnidx,respidx,ii));
              fullvpred(:,lagidx)=hinge3(nlparms,vpred(vnnidx,respidx,ii));
              
              laggain(lagidx,1:2,ii,nlidx)=nlparms;
           end
           
          case {2,3,4,5},
           % separate pos and neg gain, single threshold
           
           % fit each gain first
           for lagidx=1:fixlen,
              
              if nlidx==2,
                 % causal -- previous fixation frame
                 r1p=rpred(:,resp0idx,ii);
                 r1n=rpred(:,resp0idx+1,ii);
                 v1p=vpred(:,resp0idx,ii);
                 v1n=vpred(:,resp0idx+1,ii);
                 
              elseif nlidx==3,
                 % anti-causal -- next fixation frame
                 r1p=rpred(:,resp0idx,ii);
                 r1n=rpred(:,resp0idx-2,ii);
                 v1p=vpred(:,resp0idx,ii);
                 v1n=vpred(:,resp0idx-2,ii);
                 
              elseif nlidx==4,
                 % true positive and negative
                 r1p=posrpred(:,resp0idx,ii);
                 r1n=negrpred(:,resp0idx,ii);
                 v1p=posvpred(:,resp0idx,ii);
                 v1n=negvpred(:,resp0idx,ii);
                 
              elseif nlidx==5,
                 % true positive and negative
                 r1p=posrpred(:,resp0idx,ii);
                 r1n=negrpred(:,resp0idx,ii).*0;
                 v1p=posvpred(:,resp0idx,ii);
                 v1n=negvpred(:,resp0idx,ii).*0;
                 
              elseif nlidx==6,
                 % fake positive and negative
                 r1p=fakeposrpred(:,resp0idx,ii);
                 r1n=fakenegrpred(:,resp0idx,ii);
                 v1p=fakeposvpred(:,resp0idx,ii);
                 v1n=fakenegvpred(:,resp0idx,ii);
                 
              end
              
              r0=eresp(:,lagidx);
              lagnnidx=find(~isnan(r0) & ~isnan(r1p) & ~isnan(r1n));
              
              if sum(abs(r1p(lagnnidx)))==0 & sum(abs(r1n(lagnnidx)))==0,              
                 nlparms=[0 0 0 0];
              elseif sum(abs(r1p(lagnnidx)))==0,              
                 nlparms=fithinge3(r1n(lagnnidx),r0(lagnnidx),0,1);
                 nlparms=[0 0 nlparms];
              elseif sum(abs(r1n(lagnnidx)))==0,
                 nlparms=fithinge3(r1p(lagnnidx),r0(lagnnidx),0,1);
                 nlparms=[nlparms 0 0];
              else
                 nlparms=fithinge3([r1p(lagnnidx) r1n(lagnnidx)],r0(lagnnidx),0,2);
              end
              
              fullrpred(:,lagidx)=hinge3(nlparms,[r1p(nnidx) r1n(nnidx)]);
              fullvpred(:,lagidx)=hinge3(nlparms,[v1p(vnnidx) v1n(vnnidx)]);
              
              laggain(lagidx,1:length(nlparms),ii,nlidx)=nlparms;
           end
         end
         
         if sum(abs(nlparms))>0,
            tr=eresp(nnidx,:);
            goodidx=find(~isnan(fullrpred(:)) & ~isnan(tr(:)));
            
            nlparm=fithinge3(fullrpred(goodidx),tr(goodidx));
            fullrpred(:)=hinge3(nlparm,fullrpred(:));
            fullvpred(:)=hinge3(nlparm,fullvpred(:));
         end
         
         tp=fullrpred(:);
         goodidx=find(~isnan(tp) & ~isnan(tr(:)));
         if std(tp(goodidx))>0,
            xc(1,ii,nlidx)=xcov(tp(goodidx),tr(goodidx),0,'coeff');
         end
         
         tvp=fullvpred(:);
         tvr=vresp(vnnidx,:);
         vgoodidx=find(~isnan(tvp) & ~isnan(tvr(:)));
         
         if std(tvp(:))>0,
            allpredxc(1,ii,nlidx)=xcov(tvp(vgoodidx),tvr(vgoodidx),0,'coeff');
         end
         
         allrpred(:,:,ii,nlidx)=fullrpred;
         allvpred(:,:,ii,nlidx)=fullvpred;
         
      end
      fprintf('.');
   end
   
   dbsetqueue;
   
end
fprintf('\n');

finalvpred=zeros(length(vnnidx),fixlen,nlcount);
finalrpred=zeros(length(nnidx),fixlen,nlcount);

% save results
for nlidx=1:nlcount,
   
   g=[0.1 0.8 0.1];
   
   % force all nlidxs>3 to have same sfsidx
   % force all nlidxs to have sfsidx for nlidx=3
   if 1,
      sfsmaxidx=showidx1;
   elseif 1 | nlidx>1,
      txc=conv2(xc(1,:,3),g,'same');
      sfsmaxidx=min(find(txc==max(txc)));
   else
      txc=conv2(xc(1,:,nlidx),g,'same');
      sfsmaxidx=min(find(txc==max(txc)));
   end
   
   sfsrec(nlidx)=sfsmaxidx;
   
   fprintf('resp %d:  sfsmaxidx=%d\n',nlidx,sfsmaxidx);
   
   strf(nlidx).h=thf(:,:,sfsmaxidx);
   strf(nlidx).nlparms=[0];
   strf(nlidx).nltype='';
   strf(nlidx).mS=mSall;
   strf(nlidx).parms.kernfmt=params.kernfmt;
   strf(nlidx).parms.iconside=params.iconside;
   strf(nlidx).name=sprintf('%s nlidx=%d, sfs=%d',params.cellid,...
                            nlidx,sfsmaxidx);
   
   tvp=allvpred(:,:,sfsmaxidx,nlidx);
   tvp=tvp(:);
   tvr=vresp(vnnidx,:);
   vgoodidx=find(~isnan(tvp) & ~isnan(tvr(:)));
   [cxy,exy,tt,p]=randxcov(tvp(vgoodidx),tvr(vgoodidx),0,100);
   
   predxc(nlidx)=cxy(1);
   prederr(nlidx)=exy(1);
   predp(nlidx)=p;
   
   tvp=nansum(allvpred(:,:,sfsmaxidx,nlidx)')';
   tvr=varesp(vnnidx);
   vgoodidx=find(~isnan(tvp) & ~isnan(tvr(:)));
   [cxy,exy,tt,p]=randxcov(tvp(vgoodidx),tvr(vgoodidx),0,100);
   
   predfix(nlidx)=cxy(1);
   predfixerr(nlidx)=exy(1);
   
   finalrpred(:,:,nlidx)=allrpred(:,:,sfsmaxidx,nlidx);
   finalvpred(:,:,nlidx)=allvpred(:,:,sfsmaxidx,nlidx);
end

modeloptridx=ones(nlcount,1);
modeloptsfs=sfsrec;


% clear random stuff
clear tstim stim0 cstim0 cxy exy tt p tvp tvr vgoodidx sfsmaxidx
clear tp goodidx tr nlparm pfullrpred nfullrpred pfullvpred nfullvpred
clear r0 r1 r1p r1n allrpred

return


