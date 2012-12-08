% function res=kvatune(cellid,batch,showextra)
%
% load results from kernfile in sRunData and display tuning
% info ... kernfile generated from kernvsatt or xcdms
%
% showextra - default = 'none'
%
% created SVD 10/31/04 - hacked from dmsres.m
%
function res=kvedtune(runidx,batch,showextra)

dbopen;
if isnumeric(runidx),
   sql=['SELECT * from sRunData WHERE id=',num2str(runidx)];
   rundata=mysql(sql);
   batchcount=length(rundata);
   cellid=rundata.cellid;
else
   cellid=runidx;
   goodbatch=zeros(1,length(batch));
   batchcount=0;
   for ii=1:length(batch),
      sql=['SELECT * from sRunData WHERE cellid="',cellid,'"',...
           ' AND batch=',num2str(batch(ii))];
      trd=mysql(sql);
      if ~isempty(trd),
         goodbatch(ii)=1;
         batchcount=batchcount+1;
         rundata(batchcount)=trd;
      end
   end
   batch=batch(find(goodbatch));
end

if batchcount==0,
   disp('no entries found in db!');
   if nargout>0,
      r=0;
   end
   return
end

if ~exist('showextra','var'),
   showextra='none';
end
showextra=strsep(showextra,':');
if length(showextra)>1,
   prefframe=showextra{2};
   showextra=showextra{1};
end

global GCOLORMAP
GCOLORMAP=redblue;

rcsetstrings;

kernfile=[rundata(1).respath,rundata(1).kernfile,'.gz'];
if ~exist(kernfile,'file'),
   USERESFILE=1;
   kernfile=[rundata(1).respath,rundata(1).resfile,'.gz'];
else
   USERESFILE=0;
end

if ~exist(kernfile,'file'),
   fprintf('%s not found!\n',kernfile);
   res.cellid=cellid;
   res.cc=nan;
   res.bad=1;
   return
end

fprintf('Loading: %s\n',kernfile);
zload(kernfile);

if ~isfield(params,'bootcount'),
   if isfield(params,'fitboot') & params.fitboot>0,
      params.bootcount=params.fitboot;
   else
      params.bootcount=1;
   end
elseif params.bootcount==0,
   params.bootcount=1;
end

% raw kernels
if USERESFILE,
   H=zeros([size(strf(params.nlidxsave,1).h,1),params.bootcount]);
   mS=zeros([size(strf(params.nlidxsave,1).mS,1),params.bootcount]);
   nlparms=zeros(1,params.bootcount,attcount);
   for bootidx=1:params.bootcount,
      H(:,bootidx)=sum(strf(params.nlidxsave,bootidx).h,2) ./ ...
          length(strf(params.nlidxsave,bootidx).tempresp);
      mS(:,bootidx)=strf(params.nlidxsave,bootidx).mS;
      nlparms(bootidx)=-strf(params.nlidxsave,bootidx).nlparms;
   end
   mH=mean(H,2);
   for bootidx=1:params.bootcount,
      if H(:,bootidx)'*mH < 0,
         H(:,bootidx)=-H(:,bootidx);
      end
   end
   
   blen=0;
   randxc=predres(end).prederr(params.nlidxsave);
   %meanactresp=nanmean(predres(end).act_resp{1});
   
   meanactresp=nanmean(resp(:)).*72;
   choserange=params.nlidxsave:strfcount:size(predres(end).mod_psth{1},3);
   meanpredresp=nanmean(predres(end).mod_psth{1}(:,choserange));
   scale2Hz=(meanactresp./meanpredresp);
else
   H=zeros([size(vstrf(1).h,1),params.bootcount,attcount]);
   mS=zeros([size(vstrf(1).mS,1),params.bootcount,attcount]);
   nlparms=zeros(2,params.bootcount,attcount);
   for attidx=1:size(vstrf,1),
      H(:,:,attidx)=cat(2,vstrf(attidx,1,1,:).h);
      mS(:,:,attidx)=cat(2,vstrf(attidx,1,1,:).mS);
      nlparms(:,:,attidx)=cat(2,vstrf(attidx,1,1,:).nlparms);
   end
   
   meanactresp=mean(nlparms(1,:,1)).*60;
   scale2Hz=repmat(60,1,params.bootcount);
end

% if extra dc dimension was specified in pfft, set those coeffs to zero
kernfmt=params.kernfmt;
if strcmp(kernfmt,'fft') | strcmp(kernfmt,'lin2'),
   phasecount=4;
elseif length(kernfmt)>4 & ...
      (kernfmt(end-1)=='+' | kernfmt(end-1)=='-'),
   phasecount=str2num(kernfmt(end));
else
   phasecount=1;
end

chancount=size(H,1)./phasecount;

showSFOR=1;
if showSFOR & ~strcmp(params.kernfmt,'space'),
   Xmax=sqrt(chancount.*2);
   kernfmt='pfftgr';
   [cfilt,cfiltconj]=gencfilt(Xmax,Xmax);
else
   Xmax=params.stimloadparms{3};
   kernfmt=params.kernfmt;
end

SHOWALLKERNS=0;
if SHOWALLKERNS,
   % show jackknifed kernels
   figure(1);
   showkern(H,kernfmt);
   
   kshowcount=min([params.bootcount 12]);
   for ii=1:kshowcount,
      for attidx=1:attcount,
         subplot(attcount,kshowcount,(attidx-1).*kshowcount+ii),
         
         if ii==1,
            title(sprintf('%s %d/%d',cellid,...
                          vstrf(attidx,1,1,ii).parms.sfsfit,...
                          vstrf(attidx,1,1,ii).parms.sigfit));
         else
            title(sprintf('%d/%d',vstrf(attidx,1,1,ii).parms.sfsfit,...
                          vstrf(attidx,1,1,ii).parms.sigfit));
         end
      end
   end
   drawnow
end

%
% FIND BEST NCART MATCH TO KERNEL
%
% load ncart stimuli

nstimloadparms=params.stimloadparms;
nstimfilterparms=params.stimfilterparms;

% don't crop out ed stim
nstimloadparms{1}=nstimloadparms{3};
edmov=loadimfile('/auto/k5/david/sample/stim/edstim.imsm',...
                 0,0,nstimloadparms{:});
edmovlen=size(edmov,3);
%disp('only using contour stim in edmov');
%edmov=edmov(:,:,1:edmovlen/3);
%edmovlen=size(edmov,3);

% figure out mean and variance in space domain for fit
% data. normalize grating and natural image sets to match.
mm=0;
ss=0;
flentot=0;
pteststim=zeros(size(edmov(:,:,1)));

for fidx=1:filecount,
   
   nstimloadparms=params.stimloadparms;
   % deal with online stim aperture sizing
   if (strcmp(params.stimloadcmd,'loadimfile') |...
       strcmp(params.stimloadcmd,'loadimscaled')) & ...
         length(nstimloadparms)>0 & nstimloadparms{1}==0,
      nstimloadparms{1}=round(nstimloadparms{3} .* ...
                              params.stimcrfs(fidx)./params.stimwindowcrf);
   end
   
   flen=imfileinfo(params.stimfiles{fidx},1);
   if flen>(4000./filecount) & ismember(params.batch,[77,78,89]),
      flen=round(4000./filecount);
   elseif flen>20000./filecount,
      flen=round(20000./filecount);
   end
   
   teststim=feval(params.stimloadcmd,params.stimfiles{fidx},...
                  1,flen,nstimloadparms{:});
   
   teststim=reshape(teststim,Xmax*Xmax,size(teststim,3));
   mm=mm+mean(teststim(:)).*params.resplen(fidx);
   mpf=mean(teststim,1);
   sstd=std(teststim);
   sspf=(mean((teststim-repmat(mpf,Xmax*Xmax,1)).^2));
   spf=sqrt(mean((teststim-repmat(mpf,Xmax*Xmax,1)).^2))./(mpf+(mpf==0));
   ss=ss+mean(sstd) .* params.resplen(fidx);
end
mm=mm./sum(params.resplen);
ss=(ss./sum(params.resplen));
pteststim=pteststim./sum(params.resplen);

fprintf('fit set stats: mm=%.1f ss=%.1f\n',mm,ss);

% currently testing several different means of normalizing natural
% image and grating sets

emm=mean(edmov(:));
ess=mean(std(edmov(:,:)));
edmov(:,:)=((edmov(:,:)-emm)./ess).*ss+mm;

if 1,
   % adjust mm and ss of carts to make them vary as much as natural images
   disp('splitting into high and low contrast groups');
   
   meanadj=[1 1   2 2];
   stdadj= [1 1.5 1.0 1.5]; 
   %meanadj=[1 1 1 1];
   %stdadj= [1 0.75 1.5 2.25]; 
   %meanadj=[1 1.1 1.2 1.3];
   %stdadj= [1 1 1 1];
   sc=length(meanadj);
   
   edmov=repmat(edmov,[1 1 1 sc]);
   
   for ii=2:length(meanadj),
      edmov(:,:,:,ii)=(edmov(:,:,:,ii)-emm).*stdadj(ii)+emm*meanadj(ii);
   end
   edmovlen=edmovlen*4;
   edmov=reshape(edmov,Xmax,Xmax,edmovlen);
   emm=mean(edmov(:));
   ess=mean(std(edmov(:,:)));
   
   % re-normalize
   edmov(:,:)=((edmov(:,:)-emm)./ess).*ss+mm;
end

edmov(edmov<0)=0; edmov(edmov>255)=255;
fprintf('final  edm: emm=%.1f ess=%.1f\n',...
        mean(edmov(:)),mean(std(edmov(:,:))));

if 0,
   
   hist(mean(edmov(:,:))) 
   hist(std(edmov(:,:))) 
end

% convert normalized test stimuli into linearized domain for doing
% predictions with the SRF

% batch 82,99 - bg=46
% batch 77 - bg=20

if ~isempty(params.stimfiltercmd),
   edpfft=feval(params.stimfiltercmd,edmov,nstimfilterparms{:});
end

% subtract mean from pfft stim, since that's always necessary for
% any pred/corr thing:
edpfft=edpfft-repmat(mS(:,1),[1 edmovlen]);

cct=size(H,1);
Hnorm=H(:,:,1);
Hnorm=(Hnorm-repmat(mean(Hnorm),[cct 1]));

HNcount=size(Hnorm,2);

edmatchcc=zeros(edmovlen,HNcount);

fprintf('computing ccs');
for bootidx=1:HNcount,
   fprintf('.');
   edmatchcc(:,bootidx)=...
       (edpfft-repmat(mean(edpfft),[cct 1]))'*Hnorm(:,bootidx) ./ ...
       sqrt(sum((edpfft-repmat(mean(edpfft),[cct 1])).^2)');
end
fprintf('\n');


% predict response to each natural image
edmatch=zeros(edmovlen,params.bootcount);

fprintf('testing edmov resp...\n');
for bootidx=1:params.bootcount,
   edmatch(:,bootidx)=edpfft'*H(:,bootidx,1);
end
edmatch=edmatch+repmat(nlparms(1,:,1),[edmovlen 1]);

% this is an old, unsuccessful attempt at more realstic mean firing
% rate measurements.
if USERESFILE,
   tr=[edmatch];
   scale2Hz=repmat(72,[1 size(tr,2)]); % meanactresp./mean(tr);
   clear tr
end

% adjust spike rate to something resembling actual rates. hopefully
% doing so in a way that doesn't create big variance between
% jackknifed strfs
for bootidx=1:params.bootcount,
   edmatch(:,bootidx)=edmatch(:,bootidx).*scale2Hz(bootidx);
end

% figure out mean and stderr on each of the predicted responses

edmatchm=mean(edmatch,1);
edmatcherr=std(edmatch-repmat(edmatchm,[size(edmatch,1) 1]),...
                1,2).*sqrt(params.bootcount-1);
edmatch=mean(edmatch,2);

edmatchccm=mean(edmatchcc,1);
edmatchccerr=std(edmatchcc-repmat(edmatchccm,[size(edmatchcc,1) 1]),...
                1,2).*sqrt(params.bootcount-1);
edmatchcc=median(edmatchcc,2);

%
% figure out peak responses to each stim class

bootidx=1;

emin=min(find(edmatch(:,bootidx)==min(edmatch(:,bootidx))));
emax=min(find(edmatch(:,bootidx)==max(edmatch(:,bootidx))));

[dummy,edset]=sort(-edmatch(:,bootidx));
edsetlen=30;
edset=edset(1:edsetlen);
ed1=edset(1);

edocount=12;
edsepcount=11;
edpolcount=3;
[em1,em2,em3]=ind2sub([edocount edsepcount edpolcount*4],ed1);


edmax=sortrows([edmatch(:) edmatcherr(:)]);
edmaxerr=edmax(end,2);
edmax=edmax(end,1);

edccmax=sortrows([edmatchcc(:) edmatchccerr(:)]);
edccmaxerr=edccmax(end,2);
edccmax=edccmax(end,1);

edcc1=find(edmatchcc(:)==edccmax);
[ecc1,ecc2,ecc3]=ind2sub([edocount edsepcount edpolcount*4],edcc1);


fprintf('edmov max: %5.1f +/- %.1f  cc: %.2f +/- %.2f\n',...
        edmax,edmaxerr,edccmax,edccmaxerr);

titles{1}=sprintf('%s: m %.1f Hz',cellid,meanactresp);
titles{2}=sprintf('ed: %.1f +%.1f (%d,%d,%d)',edmax,edmaxerr,em1,em2,em3);

% fourier domain optimal gratings
feopt=edpfft(:,edset(1));
feset=edpfft(:,edset);

%mean(H(:,:,1),2)
tk=cat(3,mean(H(:,:,1),2),feopt);

eopt=edmov(:,:,edset(1))./255;

% display optimal grating from each class and optimal natural image
figure(1);
clf
showkern(repmat(tk,[1,2]),kernfmt,[Xmax Xmax],titles,1);

subplot(2,2,3);
cla
imagesc(repmat(eopt,[1 1 3])); axis image; axis off;
title(titles{2});

fullpage portrait

drawnow

if ~strcmp(showextra,'none'),
   figure(2);
   clf
end

if strcmp(showextra,'curves'),
   % show marginal tuning curves
   disp('computing tuning curves');
   obincount=16;
   sfbincount=8;
   ortuning=zeros(obincount,params.bootcount);
   sftuning=zeros(sfbincount,params.bootcount);
   orslice=zeros(obincount,params.bootcount);
   sfslice=zeros(sfbincount,params.bootcount);
   seprat=zeros(params.bootcount,1);
   cm=zeros(params.bootcount,1);
   cmstd=zeros(params.bootcount,1);
   
   for bootidx=1:params.bootcount
      tsfIR=pfft2sf(H(:,bootidx,1),params.kernfmt);
      [tsfgr,obins,sfbins]=sf2gr(tsfIR,obincount,sfbincount);
      psfgr=tsfgr.*(tsfgr>0);
      
      if bootidx==1,
         tsf0=pfft2sf(mean(H(:,:,1),2));
         tsf0=sf2gr(tsf0,obincount,sfbincount);
         mor=mean(tsf0,2);
      end
      
      [u,s,v]=svd(tsfgr);
      if sum(diag(s))>0,
         seprat(bootidx)=s(1)./sum(diag(s));
      end
      if 0,
         if u(:,1)'*mor >0,
            por=u(:,1).*s(1);
            psf=v(:,1).*s(1);
         else
            por=-u(:,1).*s(1);
            psf=-v(:,1).*s(1);
         end
         mor=mor+por;
      else
         por=mean(tsfgr,2);
         psf=mean(tsfgr',2);
      end
      
      % adjust circstats to deal with extra pi half of circle
      [cm(bootidx),cstd(bootidx)]=circstats(por(:,1));
      cstd(bootidx)=cstd(bootidx)./2 .*180/pi; 
      cm(bootidx)=cm(bootidx); % ./2 .*180/pi;
      
      ortuning(:,bootidx)=por(:,1);
      sftuning(:,bootidx)=psf(:,1);
      
      peaksfidx=min(find(psf(:,1)==max(psf(:,1))));
      orslice(:,bootidx)=tsfgr(:,peaksfidx);
      peakoridx=min(find(por(:,1)==max(por(:,1))));
      sfslice(:,bootidx)=tsfgr(peakoridx,:)';
   end
   [mcm,ecm]=circstats(ones(size(cm)),cm);
   mcm=mod(mcm.*180/pi/2,180);
   ecm=ecm.*180/pi/2.*sqrt(params.bootcount-1);
   mcstd=mean(cstd);
   ecstd=std(cstd,1).*sqrt(params.bootcount-1);
   mseprat=mean(seprat);
   eseprat=std(seprat,1).*sqrt(params.bootcount-1);
   
   mor=mean(ortuning,2);
   eor=std(ortuning-repmat(mean(ortuning,1),[size(ortuning,1) 1]),1,2).*sqrt(params.bootcount-1);
   msf=mean(sftuning,2);
   esf=std(sftuning-repmat(mean(sftuning,1),[size(sftuning,1) 1]),1,2).*sqrt(params.bootcount-1);
   
   peaksfidx=find(msf==max(msf));
   morslice=tsfgr(:,peaksfidx);
   
   tsf=pfft2sf(mean(H(:,:,1),2),kernfmt);
   tsf=flipud(sf2gr(tsf,obincount,sfbincount)');
   tsfpred=flipud((mor*msf')');
   ttsf=cat(3,[tsf(:) tsf(:)],[tsfpred(:) tsfpred(:)]);
   
   showkern(ttsf,'space',[length(sfbins) length(obins)],{},1);
   colormap(redblue);
   
   subplot(2,2,1);
   title(sprintf('%s predcc=%.2f',cellid,predxc(1)));
   
   subplot(2,2,4);
   title(sprintf('sepidx=%.2f',mseprat));
   
   subplot(2,2,2);
   %errorbar(sfbins,msf,esf,'k-');
   plot(sfbins,msf,'k-','LineWidth',2);
   %plot(sfbins,mean(sfslice,2),'k-','LineWidth',2);
   hold on
   plot([sfbins(1) sfbins(1)],max(msf)+mean(esf).*[1 -1],'k-','LineWidth',2);
   sfedge=(sfbins(2)-sfbins(1))./2;
   plot([sfbins(1)-sfedge sfbins(end)+sfedge],[0 0],'k--');
   hold off
   set(gca,'YTickLabel',[]);
   axis([sfbins(1)-sfedge sfbins(end)+sfedge ...
         min([msf-mean(esf); -mean(esf)])*1.1 ...
         max([msf+mean(esf); mean(esf)]).*1.1]);
   xlabel('spatial freq (cyc/RF)');
   ylabel('rel response');
   
   subplot(2,2,3);
   plot(obins,mor,'k-','LineWidth',2);
   %plot(obins,mean(orslice,2),'k-','LineWidth',2);
   hold on
   plot([obins(1) obins(1)],max(mor)+mean(eor).*[1 -1],'k-','LineWidth',2);
   oedge=(obins(2)-obins(1))./2;
   plot([obins(1)-oedge obins(end)+oedge],[0 0],'k--','LineWidth',1);
   hold off
   set(gca,'YTickLabel',[]);
   axis([obins(1)-oedge obins(end)+oedge ...
         min([mor-mean(eor); -mean(eor)])*1.1 ...
         max([mor+mean(eor); mean(eor)]).*1.1]);
   title(sprintf('cm=%.1f sem=%.1f cstd=%.1f',mcm,ecm,mcstd));
   xlabel('orientation (deg)');
   ylabel('rel response');
   
   set(gcf,'PaperOrientation','portrait','PaperPosition',[0.75 2.5 7 4.5]);
   
elseif strcmp(showextra,'ed');
   
   concount=sc;
   edmov=reshape(edmov,Xmax,Xmax,edocount,edsepcount,edpolcount,concount);
   edmatch=reshape(edmatch,edocount,edsepcount,edpolcount,concount);
   edmatchcc=reshape(edmatchcc,edocount,edsepcount,edpolcount,concount);
   
   % convexity
   mm=squeeze(mean(mean(mean(edmatch,1),2),4));
   mmcc=squeeze(mean(mean(mean(edmatchcc,1),2),4));
   fprintf('convexity vs resp:\n');
   fprintf('convex: %.1f (%.3f)\n',mm(2),mmcc(2));
   fprintf('outline: %.1f (%.3f)\n',mm(1),mmcc(1));
   fprintf('concave: %.1f (%.3f)\n',mm(3),mmcc(3));
   
   cvset=1;
   
   % acuteness  sepidx: 30 (1,11) 90 (3,9) 150 (5,7)
   mm=squeeze(mean(mean(mean(edmatch(:,[1 3 5],cvset,:),1),3),4)+...
              mean(mean(mean(edmatch(:,[11 9 7],cvset,:),1),3),4));
   mmcc=squeeze(mean(mean(mean(edmatchcc(:,[1 3 5],cvset,:),1),3),4)+...
                mean(mean(mean(edmatchcc(:,[11 9 7],cvset,:),1),3),4));
   fprintf('acuteness vs resp:\n');
   fprintf('30: %.1f (%.3f)\n',mm(1),mmcc(1));
   fprintf('90: %.1f (%.3f)\n',mm(2),mmcc(2));
   fprintf('150: %.1f (%.3f)\n',mm(3),mmcc(3));
   
   % sharp (1:5), smooth (6:11)
   mm=[mean(mean(mean(mean(edmatch(:,1:5,cvset,:))))); ...
       mean(mean(mean(mean(edmatch(:,7:11,cvset,:)))))];
   mmcc=[mean(mean(mean(mean(edmatchcc(:,1:5,cvset,:))))); ...
         mean(mean(mean(mean(edmatchcc(:,7:11,cvset,:)))))];
   fprintf('sharp vs smooth:\n');
   fprintf('sharp: %.1f (%.3f)\n',mm(1),mmcc(1));
   fprintf('smooth: %.1f (%.3f)\n',mm(2),mmcc(2));
   
   % display best ed stim
   estim=edmov(:,:,edset);
   estim=reshape(estim,Xmax*Xmax,10,length(edset)/10);
   
   showkern(estim,'space');
   
   keyboard
end

% determine how many pixels wide the RC window is:
if params.stimloadparms{1}==0,
   stimdiamrat=params.stimwindowcrf./params.stimcrfs(1);
else
   stimdiamrat=params.stimloadparms{3}./params.stimloadparms{1};
end

sql=['SELECT * FROM sCellFile',...
     ' WHERE cellid="',cellid,'"',...
     ' AND concat(stimpath,stimfile)="',params.stimfiles{1},'"'];
cellfiledata=mysql(sql);

windowpix=stimdiamrat.*cellfiledata(1).stimwindowsize;
windowdeg=rfdiam(cellid,windowpix);
fprintf('stim window used (deg) = %.2f\n',windowdeg);

drawnow


% if nargout==0, don't return anything in res
clear res

if nargout>0,
   
   % if requested output useful results to res 
   res.cellid=cellid;
   res.bad=0;
   res.windowdeg=windowdeg;
   res.cc=cc;
   res.H=H;
   res.ccraw=ccraw;
   res.meanactresp=meanactresp;
   res.meanpredresp=mean(natmatch(:));
   res.cartmin=cartmatch(cmin);
   res.cartmax=cartmax;
   res.cartparms=[cartfrange(cm1),cartorange(cm2),cartprange(cm3)];
   res.polmin=polmatch(pmin);
   res.polmax=polmatch(pm1,pm2,pm3,1);
   res.polparms=[polprange(pm3),polcrange(pm1),polrrange(pm2)];
   res.hypmin=hypmatch(hmin);
   res.hypmax=hypmatch(hm1,hm2,hm3,1);
   res.hypparms=[hyporange(hm1),hypprange(hm3),hypfrange(hm2)];
   res.natparms=nmax;
   res.edparms=[em1 em2 em3];
   res.edccparms=[ecc1 ecc2 ecc3];
   res.pnrat=pnrat;
   res.natfracbetter=natfracbetter;
   res.natovercart=natovercart;
   res.natmin=[natmin1 natmin2 nsort(1)];
   res.natmax=[natmax1 natmax2 nsort(end)];
   res.natmaxerr=[natmax1err natmax2err natmaxerr];
   
   res.cartvs=[res.cartmax ncartmax natmax1 edmax];
   res.cartvserr=[cartmaxerr ncartmaxerr natmax1err edmaxerr];
   res.cartvscc=[cartccmax ncartccmax natccmax edccmax];
   res.cartvsccerr=[cartccmaxerr ncartccmaxerr natccmaxerr edccmaxerr];
   res.stimclassstd=[cartccmax ncartccmax natccmax edccmax];
   res.synthvs=[synthmax natmax2];
   res.synthvserr=[synthmaxerr natmax2err];
   
   res.sparseness=[sparseness(cartmatch(:)-mmmin) ...
                   sparseness(polmatch(:)-mmmin) ...
                   sparseness(hypmatch(:)-mmmin) ...
                   sparseness(natmatch(:)-mmmin)];
   res.sparseplus=[sparseness(cartmatch(:)-mmmin) ...
                   sparseness([polmatch(:);hypmatch(:)]-mmmin) ...
                   sparseness([cartmatch(:);polmatch(:);hypmatch(:)]-mmmin) ...
                   sparseness(natmatch(:)-mmmin)];
   
   res.entropy=[hcart hpol hhyp hnat];
   
   res.randxc=randxc;
   res.beta1=cat(2,allbeta{:,end-1});
   res.beta2=cat(2,allbeta{:,end});
end

%keyboard

return
