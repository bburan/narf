% function res=kvadres(cellid,batch)
%
% load results from kernfile and display attentional modulation info ... 
% for kernfile generated from kernvsattd.m
%
% created SVD 9/18/04 - hacked from kvares.m
%
function res=kvadres(runidx,batch)

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

global GCOLORMAP
GCOLORMAP=redblue;

rcsetstrings;

if length(rundata)>1,
   disp('cellres-kerncomp: only displaying first batch');
end

fprintf('kvadres.m: loading %s\n',...
        [rundata(1).respath,rundata(1).kernfile,'.gz']);
zload([rundata(1).respath,rundata(1).kernfile,'.gz']);

% show spatial kernels for each attention condition and jackknife

kernfmt=params.kernfmt;
if strcmp(kernfmt,'pfftgr'),
   kernfmt='pfft';
   iconside=[sqrt(spacecount*2) sqrt(spacecount*2)];
end

H0=zeros([size(vstrf(1).h,1),params.bootcount,attcount]);
Hd=zeros([size(vstrf(1).h,1),params.bootcount,attcount]);
mS=zeros([size(vstrf(1).mS,1),params.bootcount,attcount]);

for attidx=1:attcount,
   H0(:,:,attidx)=cat(2,vstrf(1,attidx,1,:).h);
   Hd(:,:,attidx)=cat(2,vstrf(2,attidx,1,:).h);
   mS(:,:,attidx)=cat(2,vstrf(1,attidx,1,:).mS);
end

if nargout==0,
   
   figure(1);
   showkern(H0,kernfmt,iconside,{},0);
   
   kshowcount=min([params.bootcount 12]);
   for ii=1:kshowcount,
      for attidx=1:attcount,
         subplot(attcount,kshowcount,(attidx-1).*kshowcount+ii),
         
         if attidx==1 & ii==1,
            title(sprintf('%s %d/%d',cellid,vstrf(1,attidx,1,ii).parms.sfsfit,...
                          vstrf(1,1,1,ii).parms.sigfit));
         else
            title(sprintf('%d/%d',vstrf(1,attidx,1,ii).parms.sfsfit,...
                          vstrf(1,attidx,1,ii).parms.sigfit));
         end
         if ii==1,
            ylabel(sprintf('pidx=%d',attidx));
         end
      end
      xlabel(sprintf('boot=%d',ii));
   end
else
   disp('nargout>0: skipping figures');
end

targpredsraw=zeros(size(mpatches,2),params.bootcount,attcount);
for bootidx=1:params.bootcount,
   for attidx=1:attcount
      ts=mpatches-repmat(mS(:,bootidx,attidx),[1 size(mpatches,2)]);
      
      targpredsraw(:,bootidx,attidx)=ts'*H0(:,bootidx,attidx);
      targpredsraw(:,bootidx,attidx)=...
          dcgain(vstrf(1,attidx,1,bootidx).nlparms,...
                 targpredsraw(:,bootidx,attidx));
   end
end

H0=mean(H0,2);
Hd=mean(Hd,2);
mS=mean(mS,2);

H2=repmat(H0,[1 3]);
dcgp=ones(2,attcount);
dcgp(:,1)=vstrf(1,1,1,1).nlparms;

H2(:,2,:)=H0+Hd./2;
H2(:,3,:)=H0-Hd./2;

figure(2);
showkern(H2,kernfmt);

keyboard


% how similar each target is to the strf
targpreds=(H2(:,:)'* ...
    (mpatches(:,2:end)-repmat(mS(:,1),1,size(mpatches,2)-1)))' + ...
    repmat(dcgp(1,:),[size(mpatches,2)-1 1]);

% how similar each stimulus is to the strf
strfsim=(H2(:,:)'* (bstim'-repmat(mS(:,1),1,blen)))' + ...
        repmat(dcgp(1,:),[blen 1]);


if nargout>0,
   res.kernfile=[rundata(1).respath,rundata(1).kernfile,'.gz'];
   
   % save mean response
   res.mresp=nanmean(squeeze(bresp));
   
   % save predxc
   res.predxc=predxc;
   res.pxc=pxc;
   res.predxccross=squeeze(predxccross(1,:,:));
   res.predxccrossrand=squeeze(predxccross(2,:,:));
   res.randxc=randxc;
   
   % predicted baseline response to targets, mean response, rank of targets
   res.strfresp=targpreds;
   res.strfmean=mean(strfsim);
   for aa=1:size(strfsim,2),
      tt=sort(strfsim(:,aa));
      %res.strfrank(1,aa)=max([0; find(targpreds(1,aa)>tt)])./length(tt);
      %res.strfrank(2,aa)=max([0; find(targpreds(2,aa)>tt)])./length(tt);
      res.strfrank(:,aa)=(targpreds(:,aa)-mean(tt))./std(tt);
   end
   
   return
end

ap=fpatches;
%tp=mpatches;
tp=mpatches-repmat(mS(:,1),[1 size(mpatches,2)]);
ap = ap ./ repmat(mean(ap,1),size(ap,1),1);
tp = tp ./ repmat(mean(abs(tp),1),size(tp,1),1);
mp=repmat(mean(ap,2),1,size(tp,2));
%tp=tp-mp;

if size(tp,1)>spacecount,
   targmtx=rand(spacecount,1,attcount);
elseif size(tp,2)<attcount,
   targmtx=reshape(tp,spacecount,1,size(tp,2));
   targmtx(:,1,size(tp,2)+1:attcount)=0;
   targmtx(1,1,size(tp,2)+1:attcount)=1;
   targmtx=abs(targmtx).^0.5 .* sign(targmtx);
else
   targmtx=reshape(tp(:,1:attcount),spacecount,1,attcount);
   targmtx=abs(targmtx).^0.5 .* sign(targmtx);
end

for attidx=1:attcount,
   targmtx(:,:,attidx)=targmtx(:,:,attidx) ./ ...
       max(abs(targmtx(:,:,attidx))) .* max(max(abs(H2(:,:,attidx))));
end

eigmtx=cat(2,targmtx,H,H2,H);

fprintf('pred xc : indiv att state MSE (att/rand):\n');
fprintf('attidx     :         ');
fprintf(' %11d',1:attcount);
fprintf('\n');

fprintf('framecount :         ');
fprintf(' %11d',ncount);
fprintf('\n');

for nlidx=1:nlcount,
   fprintf('%-11s:',nlnames{nlidx});
   fprintf('%5.2f/%5.2f',predxc(1:2,nlidx));
   fprintf(' (p<%4.2f):',pxc(nlidx));
   fprintf(' %5.2f/%5.2f',predxccross(1:2,nlidx,2:attcount));
   fprintf('\n');
end

%xcboot(find(isnan(xcboot)))=0;
%[mean(xcboot(1,:,:),3); mean(mean(xcboot(2:end,:,:),3))]

figure(2);
eigmtx(:,:,size(eigmtx,3)+1)=nan;
showkern(eigmtx,kernfmt,iconside,{},0);
colcount=size(eigmtx,2);
rowcount=size(eigmtx,3);

maxr=max(bresp(:,1));
xx=linspace(-maxr,maxr,100);

subplot(rowcount,colcount,1);
plot(targpreds(:,1),'k--');
hold on
plot(targpreds(:,2:end));
hold off
a=axis;
axis([0 attcount+2 a(3) a(4)]);
sleg={'a'};
for aa=2:attcount,
   sleg{aa}=sprintf('%d',aa-1);
end
legend(sleg);

for ii=1:size(eigmtx,3),
   subplot(rowcount,colcount,ii*colcount);
   cla; axis off;
end

subplot(rowcount,colcount,colcount);
globaldc=globaldc.*60;
mdc=mean(globaldc,2);
edc=std(globaldc,1,2) .* sqrt(params.resampcount-1);
plot(linspace(0,attcount,params.bootcount),globaldc','--');
hold on
errorbar(mdc,edc,'ko');
hold off
title('att vs dc');

subplot(rowcount,colcount,colcount*2);
mgain=mean(globalgain,2);
egain=std(globalgain,1,2) .* sqrt(params.resampcount-1);

plot(linspace(0,attcount,params.bootcount),globalgain','--');
hold on
errorbar(mgain,egain,'ko');
hold off
title('att vs gain');
xlabel('attidx');

subplot(rowcount,colcount,2);
title(sprintf('%s no-att kernel',cellid));
for attidx=2:attcount,
   subplot(rowcount,colcount,attidx*colcount-colcount+1);
   ylabel(sprintf('attidx=%d',attidx-1));
end


subplot(rowcount,1,rowcount);
axis off;
axis ij;

details={};
details{1}=sprintf('%s: pred xc : indiv att state MSE (att/rand):',cellid);
details{2}=[sprintf('attidx     :         '),...
            sprintf(' %11d',1:attcount)];

details{3}=[sprintf('framecount :         '),...
            sprintf(' %11d',ncount)];

for nlidx=1:nlcount,
   details{3+nlidx}=[...
      sprintf('%-11s:',nlnames{nlidx}),...
      sprintf('%5.2f/%5.2f',predxc(1:2,nlidx)),...
      sprintf(' (p<%4.2f):',pxc(nlidx)),...
      sprintf(' %5.2f/%5.2f',predxccross(1:2,nlidx,2:attcount))];
end

ht=text(0,0,details,'VerticalAlignment','top','FontName','Courier');
fullpage('portrait');

if params.runclassid==3 | params.runclassid==10,
   if params.runclassid==10,
      targlist=1:size(bigpatches,3);
   end
   
   mmin=min(min(min(bigpatches(:,:,targlist))));
   mmax=max(max(max(bigpatches(:,:,targlist))));
   
   figure(3);
   clf
   for ii=1:length(targlist),
      subplot(1,length(targlist),ii);
      imagesc(bigpatches(:,:,targlist(ii)),[mmin mmax]);
      title(sprintf('%s target %d (id %d)',cellid,ii,targlist(ii)));
      axis off
      axis image
   end
   colormap(gray);
   set(gcf,'PaperOrientation','Portrait','PaperPosition',[0.25 0.25 8 10.5]);
end

%
% compare preds by different att states
%
if exist('attxc','var'),
   disp('local att preds');
   for nlidx=1:nlusecount,
      fprintf('nlidx=%d:\n',nluse(nlidx));
      for attidx=2:attcount,
         pin=attxc(attidx,attidx,nlidx);
         pout=nanmean(attxc([2:attidx-1 attidx+1:end],attidx,nlidx));
         fprintf('attidx=%d: %6.3f/%6.3f',attidx,pin,pout);
         
         for att2=2:attcount,
            if 0 & att2==attidx,
               fprintf('  x.xxx/ x.xxx');
            else
               p1=valattpreds(:,attidx,nlidx);
               p2=valattpreds(:,att2,nlidx);
               a1idx=find(~isnan(bresp(anyokidx,1,attidx)));
               a2idx=find(~isnan(bresp(anyokidx,1,att2)));
               
               pin=xcov(bresp(anyokidx([a1idx;a2idx]),:,1),...
                        [p1(a1idx);p2(a2idx)],0,'coeff');
               pout=xcov(bresp(anyokidx([a1idx;a2idx]),:,1),...
                         [p2(a1idx);p1(a2idx)],0,'coeff');
               
               fprintf(' %6.3f/%6.3f',pin,pout);
            end
         end
         fprintf('\n');
      end   
   end
end

tstrf=vstrf(1);
clear tunedata
for attidx=1:attcount,
   tstrf.h=H2(:,:,attidx);
   tunedata(attidx)=kern2tune(tstrf);
end

figure(4);
clf
obins=tunedata(1).obins;
sfbins=tunedata(1).sfbins;

patorsf=sf2gr(mpatches,length(obins),length(sfbins),1,0,kernfmt);
for attidx=1:attcount
   
   if attidx>1,
      subplot(attcount+1,4,attidx*4-3);
      
      ta=patorsf(:,:,attidx);
      ta=ta.*(ones(length(obins),1)*sqrt(sfbins));
      ta=ta.^2;
      
      mrange=max(abs(ta(:)));
      imagesc(obins,sfbins,ta',[-mrange mrange]);
      axis xy
      ylabel(sprintf('attidx=%d or',attidx));
      if attidx==2,
         title(sprintf('%s targets',cellid));
      end
   end
   
   subplot(attcount+1,4,attidx*4-2);
   mrange=max(max(abs(tunedata(attidx).orsf(:,:,1))));
   imagesc(obins,sfbins,tunedata(attidx).orsf(:,:,1)',[-mrange mrange]);
   
   if 1,
      xx=-2:2;
      gsf=exp(-(xx./0.8).^2/2);
      gsf=gsf./sum(gsf(:))';
      gor=exp(-(xx./1.0).^2/2);
      gor=(gor./sum(gor(:)));
      
      ta=tunedata(attidx).orsf(:,:,1)';
      astd=std(ta(:));
      
      ta=conv2(ta,gsf,'same');
      ta=cconv2(ta,gor);
      tas=size(ta);
      
      slev1=1.0;
      slev2=2.0;
      
      hold on
      contour(obins,sfbins,ta,[-astd*slev2  -astd*slev2],'k-');
      contour(obins,sfbins,ta,[-astd*slev1  -astd*slev1],'k--');
      contour(obins,sfbins,ta,[ astd*slev1   astd*slev1],'k--');
      contour(obins,sfbins,ta,[ astd*slev2   astd*slev2],'k-');
      hold off
   end
   
   axis xy
   if attidx==1,
      title(sprintf('%s spatial tuning',cellid));
   end
   
   subplot(attcount+1,4,attidx*4-1);
   plot(tunedata(attidx).obins,tunedata(attidx).por(:,1));
   ostep=(tunedata(attidx).obins(2)-tunedata(attidx).obins(1))./2;
   axis([tunedata(attidx).obins(1)-ostep ...
         tunedata(attidx).obins(end)+ostep -1 1]);
   
   if attidx==1,
      title('orientation');
   end
   
   subplot(attcount+1,4,attidx*4);
   plot(tunedata(attidx).sfbins,tunedata(attidx).psf(:,1));
   sfstep=(tunedata(attidx).sfbins(2)-tunedata(attidx).sfbins(1))./2;
   axis([tunedata(attidx).sfbins(1)-sfstep ...
         tunedata(attidx).sfbins(end)+sfstep -1 1]); 
   
   if attidx==1,
      title('spatial freq');
   end
end
colormap(redblue);

subplot(attcount+1,1,attcount+1);
axis off;
axis ij;
ht=text(0,0,details,'VerticalAlignment','top','FontName','Courier');
fullpage('landscape');


return


if OUTNLMODE==4,
   
   resplen=size(bresp,1);
   nluse=[4 5];
   nlusecount=length(nluse);
   
   attpreds=zeros(resplen,attcount,nlusecount,params.bootcount);
   attxc=zeros(attcount,attcount,nlusecount);
   
   for attidx=1:attcount,
      for nlidx=1:nlusecount,
         attpdummy=zeros(resplen,1);
         for bootidx=1:params.bootcount,
            if attidx==1,
               tstrf=vstrf(1,1,1,bootidx);
            else
               tstrf=vstrf(nluse(nlidx),attidx-1,1,bootidx);
            end
            
            tstim=bstim'-repmat(tstrf.mS,[1 resplen]);
            tpred=kernpredict(tstrf.h,tstim,1,0);
            
            tnltype=tstrf.nltype;
            tnlparms=tstrf.nlparms;
            if ~isempty(tnltype) & ~strcmp(tnltype,'none'),
               attpreds(:,attidx,nlidx,bootidx)=feval(tnltype,tnlparms,tpred);
            else
               attpreds(:,attidx,nlidx,bootidx)=tpred;
            end
            
            if nlidx==2,
               attpreds(:,attidx,nlidx,bootidx)=...
                   attpreds(:,attidx,nlidx,bootidx) + ...
                   attpreds(:,attidx,1,bootidx);
            end
            
            a0idx=vidx((round((bootidx-1)/params.bootcount*vcount)+1):...
                       round(bootidx/params.bootcount*vcount));
            attpdummy(a0idx)=attpreds(a0idx,attidx,nlidx,bootidx);
         end
         
         for att2=1:attcount,
            if att2==1,
               tokidx=anyokidx;
            else
               tokidx=find(~isnan(bresp(:,att2)));
            end
            attxc(attidx,att2,nlidx)=...
                xcov(attpdummy(tokidx),bresp(tokidx,1),0,'coeff');
         end
      end
   end
   attpreds(:,2:end,2,:)=attpreds(:,2:end,1,:)+attpreds(:,2:end,2,:);
   attpreds=attpreds(anyokidx,:,:,:);
else
   
   
   resplen=length(anyokidx);
   nluse=[3 6];
   nlusecount=length(nluse);
   
   attpreds=zeros(resplen,attcount,attcount-1,nlusecount,params.bootcount);
   
   for attidx=1:attcount,
      for nlidx=1:nlusecount,
         for bootidx=1:params.bootcount,
            if attidx==1,
               tstrf=vstrf(1,1,1,bootidx);
            else
               tstrf=vstrf(nluse(nlidx),attidx-1,1,bootidx);
            end
            
            tstim=bstim(anyokidx,:)'-repmat(tstrf.mS,[1 resplen]);
            tpred=kernpredict(tstrf.h,tstim,1,0);
            
            tnltype=tstrf.nltype;
            for att2idx=1:attcount-1,
               tnlparms=tstrf.nlparmsalt{att2idx};
               if ~isempty(tnltype) & ~strcmp(tnltype,'none'),
                  attpreds(:,attidx,att2idx,nlidx,bootidx)=feval(tnltype,tnlparms,tpred);
               else
                  attpreds(:,attidx,att2idx,nlidx,bootidx)=tpred;
            end
            end
         end
      end
   end
end


figure(4);
clf

plotcount=attcount-1;
linstr={'k:','b-','r-','g-','k-'};
mp=mean(attpreds,4);
for attidx=1:plotcount,
   
   aokidx=find(~isnan(bresp(anyokidx,1,attidx+1)));
   %aokidx=1:length(anyokidx);
   
   % 1H, att dcg: mp(aokidx,attidx+1,attidx,1),...
   % 4H, att dcg: mean(mp(aokidx,[2:attidx attidx+2:attcount],attidx,2),2),...
   
   if 1,
      smp=cat(2,mp(aokidx,1),...
              mp(aokidx,attidx+1,1),...
              mp(aokidx,attidx+1,2));
   else
      smp=zeros(length(aokidx),3);
      for attidx2=[1:attidx-1 attidx+1:attcount-1],
         smp(:,1)=smp(:,1) + mp(aokidx,attidx2+1,attidx2,2);
         smp(:,2)=smp(:,2) + mp(aokidx,attidx2+1,attidx,2);
      end
      smp(:,1:2)=smp(:,1:2)./(attcount-2);
      smp(:,3)=mp(aokidx,attidx+1,attidx,2);
   end
   
   tr=bresp(anyokidx(aokidx),attidx+1);
   
   if 1,
      [oo,ss]=sort(smp(:,1));
      
      smp=smp(ss,:);
      smp(:,end)=gsmooth(smp(:,end),5);
      %smp(:,2:end)=gsmooth(smp(:,2:end),3);
      
      tr=gsmooth(tr(ss),5);
   else
      [tr,ss]=sort(tr);
      
      smp=smp(ss,:);
      smp=gsmooth(smp,3);
   end
   
   subplot(2,2,attidx);
   for ii=size(smp,2):-1:1,
      plot(smp(:,ii),linstr{ii});
      hold on
   end
   %plot(tr,'r:');
   
   plot([0,round(size(smp,1)/3)],[1 1].*mean(targpredsraw(attidx+1,:,1),2),'r:');
   
   hold off
   title(sprintf('attidx=%d',attidx));
   
   axis([0 size(smp,1) 0 max(max(max(mean(attpreds,4))))]);
end

legend('noatt','dcg','local',2);



keyboard


return







figure(4);
clf
subplot(2,1,1);
plotstyle={'k:','r-','g--','b-','m--'};
attstate={'all','1','2','3','4'};
xx=linspace(nanmin(rprec(:,1)),nanmax(rprec(:,1)));
for attidx=1:attcount,
    plot(xx,sigmoid(sigparms(:,attidx),xx),plotstyle{attidx});
    hold on
end

subplot(2,1,2);
for attidx=1:attcount,
    errorbar((1:4)',sigparms(:,attidx),esigparms(:,attidx),plotstyle{attidx});
    hold on
end




if 0 & ~exist('asigparms','var'),
   spcount=3;
   asigparms=zeros(spcount,params.bootcount,attcount-1,spcount+1);
   valpred=zeros(size(bresp,1),spcount+1);
   valxc=zeros(1,size(asigparms,1)+1);
   
   % find baseline (all att) sigmoid parms
   fprintf('fitting hinge4, attidx=1');
   attidx=1;
   tokidx=find(sum(~isnan(bresp(:,1,2:end)),3)>0);
   for jj=1:params.bootcount,
      % use same jackknife sets as for fitting 
      useidx=vidx([1:round((jj-1)/params.bootcount*vcount) ...
                   round(jj/params.bootcount*vcount+1):end]);
      useidx=useidx(find(~isnan(bresp(useidx,1,attidx))));
      predidx=vidx(round((jj-1)/params.bootcount*vcount+1): ...
                   round(jj/params.bootcount*vcount));
      predidx=predidx(find(~isnan(bresp(predidx,1,attidx))));
      
      ract=bresp(useidx,1,1);
      rpred=rprec(useidx,1);
      
      tfitparms=fithinge4(rpred,ract,0);
      asigparms(:,jj,:,1)=repmat(tfitparms,[1 1 attcount-1]);
      
      valpred(predidx,1)=hinge4(asigparms(:,jj,attidx),...
                                rprec(predidx,1));
   end
   
   for sigidx=1:spcount,
      
      filler=ones(size(bresp,1),spcount);
      for fitidx=1:attcount-1,
         filler(find(~isnan(bresp(:,fitidx+1))),sigidx)=fitidx;
      end
      
      for jj=1:params.bootcount,
         % use same jackknife sets as for fitting 
         useidx=vidx([1:round((jj-1)/params.bootcount*vcount) ...
                      round(jj/params.bootcount*vcount+1):end]);
         predidx=vidx(round((jj-1)/params.bootcount*vcount+1): ...
                      round(jj/params.bootcount*vcount));
         
         ract=bresp(useidx,1,1);
         rpred=rprec(useidx,1);
         
         tfitparms=fithinge4(rpred,ract,0,filler(useidx,:));
         
         for attidx=1:attcount-1,
            useidx=[1:sigidx-1 sigidx+attidx-1 ...
                    (sigidx+attcount-1):(spcount+attcount-2)]';
            
            asigparms(:,jj,attidx,sigidx+1)=tfitparms(useidx);
            
            apredidx=predidx(find(~isnan(bresp(predidx,1,attidx+1))));
            valpred(apredidx,sigidx+1)=...
                hinge4(asigparms(:,jj,attidx,sigidx+1),rprec(apredidx,1));
         end
      end
   end
   for sigidx=1:spcount+1,
      valxc(sigidx)=xcov(bresp(anyokidx,1),valpred(anyokidx,sigidx),...
                         0,'coeff');
   end
   
   sigparms=squeeze(mean(asigparms,2));
   esigparms=squeeze(std(asigparms,1,2)).*sqrt(params.bootcount-1);
   nsigparms=sigparms./repmat(sigparms(:,1),[1 attcount-1 spcount+1]);
   nesigparms=esigparms./repmat(sigparms(:,1),[1 attcount-1 spcount+1]); 
elseif ~exist('asigparms','var'),
   spcount=4;
   asigparms=zeros(spcount,params.bootcount,attcount);
   valpred=zeros(size(bresp,1),spcount+1);
   valxc=zeros(1,size(asigparms,1)+1);
   
   % find baseline (all att) sigmoid parms
   fprintf('fitting sigmoid, attidx=1');
   attidx=1;
   tokidx=find(sum(~isnan(bresp(:,1,2:end)),3)>0);
   for jj=1:params.bootcount,
      % use same jackknife sets as for fitting 
      useidx=vidx([1:round((jj-1)/params.bootcount*vcount) ...
                   round(jj/params.bootcount*vcount+1):end]);
      useidx=useidx(find(~isnan(bresp(useidx,1,attidx))));
      predidx=vidx(round((jj-1)/params.bootcount*vcount+1): ...
                   round(jj/params.bootcount*vcount));
      predidx=predidx(find(~isnan(bresp(predidx,1,attidx))));
      
      ract=bresp(useidx,1,1);
      rpred=rprec(useidx,1);
      
      asigparms(:,jj,attidx)=fitsigmoid(rpred,ract,0);
      
      valpred(predidx,1)=sigmoid(asigparms(:,jj,attidx),...
                                 rprec(predidx,1));
   end
   
   for attidx=2:attcount,
      
      fprintf(' %d',attidx);
      
      tokidx=find(~isnan(bresp(:,1,attidx)));
      
      for jj=1:params.bootcount,
         for sigidx=1:spcount,
            % use same jackknife sets as for fitting 
            useidx=vidx([1:round((jj-1)/params.bootcount*vcount) ...
                         round(jj/params.bootcount*vcount+1):end]);
            useidx=useidx(find(~isnan(bresp(useidx,1,attidx))));
            
            predidx=vidx(round((jj-1)/params.bootcount*vcount+1): ...
                         round(jj/params.bootcount*vcount));
            predidx=predidx(find(~isnan(bresp(predidx,1,attidx))));
            
            ract=bresp(useidx,1,1);
            rpred=rprec(useidx,1);
            
            forcevalues=asigparms(:,jj,1);
            forcevalues(sigidx)=nan;
            
            tsigparms=fitsigmoid(rpred,ract,0,forcevalues);
            asigparms(sigidx,jj,attidx)=tsigparms(sigidx);
            
            valpred(predidx,sigidx+1)=sigmoid(tsigparms,rprec(predidx,1));
         end
      end
   end
   
   for sigidx=1:spcount+1,
      valxc(sigidx)=xcov(bresp(anyokidx,1),valpred(anyokidx,sigidx),...
                         0,'coeff');
   end
   
   sigparms=squeeze(mean(asigparms,2));
   esigparms=squeeze(std(asigparms,1,2)).*sqrt(params.bootcount-1);
   nsigparms=sigparms./repmat(sigparms(:,1),[1 attcount]);
   nesigparms=esigparms./repmat(sigparms(:,1),[1 attcount]);
else
   spcount=size(sigparms,1);
   
end
fprintf('\n');

figure(4);
clf

xmin=min(rprec(anyokidx,1))-std(rprec(anyokidx,1));
xmax=max(rprec(anyokidx,1))+std(rprec(anyokidx,1));
ymin=nanmin(bresp(:,1));
ymax=nanmax(bresp(:,1));
xrange=linspace(xmin,xmax,100);
plotstyle={'k:','r-','g--','b-','m--'};
attstate={'all','1','2','3','4'};
spstring={'x10','sigma','baseline','amp'};

for sigidx=1:spcount,
   subplot(spcount,1,sigidx);
   for attidx=1:attcount,
      sparms=sigparms(:,1);
      sparms(sigidx)=sigparms(sigidx,attidx);
      plot(xrange,sigmoid(sparms,xrange),plotstyle{attidx});
      hold on
   end
   
   for attidx=1:attcount,
      if attidx==1,
         tokidx=find(sum(~isnan(bresp(:,1,2:end)),3)>0);
      else
         tokidx=find(~isnan(bresp(:,1,attidx)));
      end
      
      plotcnt=min([100 length(tokidx)]);
      plot(rprec(tokidx(1:plotcnt),1),bresp(tokidx(1:plotcnt),1,1),...
           [plotstyle{attidx}(1),'.']);
   end
   hold off
   
   title(sprintf('%s att %s: valxc: %.3f->%.3f',cellid,...
                 spstring{sigidx},valxc(1),valxc(sigidx+1)));
end
legend(attstate{1:attcount});
fullpage('portrait');

if 0
subplot(2,1,2);

for attidx=2:attcount,
   
   errorbar(nsigparms(:,attidx),nesigparms(:,attidx),plotstyle{attidx});
   hold on
   
end
hold off
end


if 0
figure(4);
clf
for ii=1:attcount,
   subplot(attcount,1,ii);
   mm=max(max(abs(lth(:,:,ii))));
   if mm>0,
      imagesc(lth(:,:,ii)',[-mm mm]);
      title(sprintf('att %d local kernels',ii));
      axis image
   end
   axis off
end
colormap(gray);
fullpage('Portrait');
end


