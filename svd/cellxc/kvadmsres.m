% function z=kvadmsres(cellid,batch)
%
% load results from kernfile in sRunData and display attentional
% modulation info ... kernfile generated from kerncomp4
%
% z=data from kernfile
%
% created SVD 10/18/02 - hacked from kerncomp4res.m
%
function z=kvadmsres(runidx,batch)


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

fprintf('Loading: %s\n',[rundata(1).respath,rundata(1).kernfile,'.gz']);
zload([rundata(1).respath,rundata(1).kernfile,'.gz']);

% show spatial kernels for each attention condition and jackknife

kernfmt=params.kernfmt;
if strcmp(kernfmt,'pfftgr'),
   kernfmt='pfft';
   iconside=[sqrt(spacecount*2) sqrt(spacecount*2)];
end

H=zeros([size(vstrf(1).h,1),params.bootcount,attcount]);
H(:,:,1)=cat(2,vstrf(1,1,1,:).h);
for attidx=2:attcount,
   H(:,:,attidx)=cat(2,vstrf(end,attidx-1,1,:).h);
end

figure(1);
showkern(H,kernfmt,iconside,{},0);

kshowcount=min([params.bootcount 12]);
for ii=1:kshowcount,
   subplot(attcount,kshowcount,ii),
   title(sprintf('%d/%d',vstrf(1,1,1,ii).parms.sfsfit,vstrf(1,1,1,ii).parms.sigfit));
   
   for attidx=2:attcount,
      subplot(attcount,kshowcount,(attidx-1).*kshowcount+ii),
      title(sprintf('%d/%d',vstrf(6,attidx-1,1,ii).parms.sfsfit,...
                    vstrf(6,attidx-1,1,ii).parms.sigfit));
   end
end

H=mean(H,2);

ap=fpatches;
tp=mpatches;
ap = ap ./ repmat(mean(ap,1),size(ap,1),1);
tp = tp ./ repmat(mean(tp,1),size(tp,1),1);
mp=repmat(mean(ap,2),1,size(tp,2));
tp=tp-mp;

if size(tp,1)>spacecount,
   targmtx=rand(spacecount,1,attcount);
elseif size(tp,2)<attcount,
   targmtx=reshape(tp,spacecount,1,size(tp,2));
   targmtx(:,1,size(tp,2)+1:attcount)=0;
   targmtx(1,1,size(tp,2)+1:attcount)=1;
   targmtx=sqrt(abs(targmtx)) .* sign(targmtx);
else
   targmtx=reshape(tp(:,1:attcount),spacecount,1,attcount);
   targmtx=sqrt(abs(targmtx)) .* sign(targmtx);
end

mm=max(abs([H(:)]));  % ; globmtx(:)
for attidx=1:attcount,
   targmtx(:,:,attidx)=targmtx(:,:,attidx) ./ ...
       max(abs(targmtx(:,:,attidx))) .* max(max(abs(H(:,:,attidx))));
end
eigmtx=cat(2,targmtx,H,H);

figure(2);
showkern(eigmtx,kernfmt,iconside,{},0);
colcount=size(eigmtx,2);

nlshow=[1:4 6];
maxr=max(bresp(:,1));
xx=linspace(-maxr,maxr,100);

subplot(attcount,colcount,1);
cla; axis off;
for ii=1:size(eigmtx,3), % length(nlshow);
   subplot(attcount,colcount,ii*colcount);
   cla; axis off;
end

subplot(attcount,colcount,colcount);
globaldc=globaldc.*60;
mdc=mean(globaldc,2);
edc=std(globaldc,1,2) .* sqrt(params.resampcount-1);
plot(linspace(0,attcount,params.bootcount),globaldc','--');
hold on
errorbar(mdc,edc);
hold off
title('att vs dc');

subplot(attcount,colcount,colcount*2);
mgain=mean(globalgain,2);
egain=std(globalgain,1,2) .* sqrt(params.resampcount-1);

plot(linspace(0,attcount,params.bootcount),globalgain','--');
hold on
errorbar(mgain,egain);
hold off
title('att vs gain');
xlabel('attidx');

subplot(attcount,colcount,colcount*3);
agg=cat(1,cat(2,vstrf(6,1,1,:).nlparms),cat(2,vstrf(6,2,1,:).nlparms));
plot(agg([2 4],:)');
%plot(agg([1 3],:)');
title('local att gain');

subplot(attcount,colcount,2);
title(sprintf('%s no-att kernel',cellid));
for attidx=2:attcount,
   subplot(attcount,colcount,attidx*colcount-colcount+1);
   ylabel(sprintf('attidx=%d',attidx-1));
end

set(gcf,'PaperOrientation','Portrait','PaperPosition',[0.25 0.25 8 10.5]);

if params.runclassid==3 | params.runclassid==10,
   if params.runclassid==10,
      targlist=1:size(bigpatches,3);
   end
   
   figure(3);
   clf
   for ii=1:length(targlist),
      subplot(length(targlist),1,ii);
      imagesc(bigpatches(:,:,targlist(ii)),[0 255]);
      title(sprintf('target %d (id %d)',ii,targlist(ii)));
      axis off
      axis image
   end
   colormap(gray);
   set(gcf,'PaperOrientation','Portrait','PaperPosition',[0.25 0.25 8 10.5]);
end

lth=zeros(spacecount,params.bootcount,attcount);
lth(:,:,1)=cat(2,vstrf(1,1,1,:).h);
for ii=1:attcount-1,
   lth(:,:,ii+1)=cat(2,vstrf(end,ii,1,:).h);
   lgain=cat(2,vstrf(end,ii,1,:).nlparms);
   
   lth(:,:,ii+1)=lth(:,:,ii+1).*repmat(lgain(2,:),spacecount,1).*...
       repmat(globalgain(ii,:),spacecount,1);
end

if 1
   
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
   fprintf(' %5.2f/%5.2f',predxccross(1:2,nlidx,1:attcount-1));
   fprintf('\n');
end

%xcboot(find(isnan(xcboot)))=0;
%[mean(xcboot(1,:,:),3); mean(mean(xcboot(2:end,:,:),3))]

end





keyboard

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
else
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




%keyboard



return

% draw the strf one eigenvector at a time, and their sum
ecount=max(predspace); % min([40 spacelims]);
fprintf('ecount=%d\n',ecount);
erange=spaceuses(1:ecount);
eigtargets=ua(:,spaceuses(1:ecount))' * ...
    (mpatches-repmat(bsmean',[1 size(mpatches,2)]));

if isfield(z,'DAMPFACTOR'),
   DAMPFACTOR=DAMPFACTOR;
else
   DAMPFACTOR=1.0;
end
noiseidx=1;

Hm=squeeze(mean(lH(:,respidx,:,noiseidx,:),5));
Hmed=squeeze(median(lH(:,respidx,:,noiseidx,:),5));
Hs=squeeze(std(lH(:,respidx,:,noiseidx,:),0,5).*sqrt(bootcount));
%Hdamp=(abs(Hm)./(Hs.*DAMPFACTOR)).^2;
%Hdamp(find(Hdamp>1))=1;
Hdamp=(abs(Hm)./(Hs.*DAMPFACTOR));
Hdamp=(1-Hdamp.^(-2));
Hdamp=Hdamp.*(Hdamp>0);
Hdamp(find(isnan(Hdamp)))=0;

if 0,
   % damp out according to shrinkage filter
   H=Hm.*repmat(Hdamp(:,1),[1 attcount]);
elseif params.dodamping
   % damp out according to shrinkage filter
   fprintf('sfsidx max=%d\n',ddmax);
   H=Hm.*Hdamp;
else
   H=Hm;
end

a1mean=mean(a1(respidx,:,1,:,1),4);
a1std=std(a1(respidx,:,1,:,1),0,4).*sqrt(bootcount);
%rateconv=respfilterparms{2}(respidx)-respfilterparms{1}(respidx)
% scale c1 by 1000 to convert to Hz from spikes/ms
c1mean=mean(c1(respidx,:,1,:,1),4)*1000;
c1std=std(c1(respidx,:,1,:,1),0,4).*sqrt(bootcount)*1000;

kernmtx=ua(:,erange)*H(erange,:);
globmtx=kernmtx(:,1)*a1mean;
targmtx=reshape(mpatches(:,1:attcount),spacecount,1,attcount);
if strcmp(params.kernfmt,'pfft') | strcmp(params.kernfmt,'pfftgr'),
   targmtx(37,:)=0;  % zero out the dc term for 16x16 poower
end
mm=max(abs([kernmtx(:)]));  % ; globmtx(:)
globmtx=globmtx./max(abs(globmtx(:))).*mm;
for attidx=1:attcount,
   targmtx(:,:,attidx)=targmtx(:,:,attidx) ./ ...
       max(abs(targmtx(:,:,attidx))) .* max(max(abs(kernmtx(:,:))));

   if 1 & max(abs(kernmtx(:,attidx)))>0,
      kernmtx(:,attidx)=kernmtx(:,attidx) ./ ...
          max(abs(kernmtx(:,attidx))) .* max(abs(globmtx(:,attidx)));
   end
end
eigmtx=cat(2,targmtx,reshape(globmtx,spacecount,1,attcount),...
           reshape(kernmtx,spacecount,1,attcount));

attcc=zeros(attcount);
for xx=1:attcount,
   if sum(abs(kernmtx(:,xx)))>0,
      for yy=1:attcount,
         attcc(xx,yy)=xcov(kernmtx(:,xx),targmtx(:,yy),0,'coeff');
      end
   end
end

figure(1);
clf
kernfmt=params.kernfmt;
if strcmp(kernfmt,'pfftgr'),
   kernfmt='pfft';
   iconside=[sqrt(spacecount*2) sqrt(spacecount*2)];
end
showkern(eigmtx,kernfmt,iconside,{},0);

subplot(attcount,3,1);
title(sprintf('%s target',cellid));
subplot(attcount,3,2);
title(sprintf('r=%d global',respidx));
subplot(attcount,3,3);
title(sprintf('r=%d local',respidx));

for attidx=1:attcount,
   subplot(attcount,3,attidx*3-2);
   ylabel(sprintf('attidx=%d',attidx-1));
end
set(gcf,'PaperOrientation','Portrait','PaperPosition',[0.25 0.25 8 10.5]);

figure(2);
clf

%ecount=20;
rH=squeeze(lH(1:ecount,:,2:end,1,:));
H=mean(rH,3);
eH=std(rH,1,3) .* sqrt((resampcount-1)/resampcount);
mH=mean(H,2);

% subtract relevant means
rH=rH-repmat(H,[1 1 bootcount]);
H=H-repmat(mH,[1 (attcount-1)]);

Sb=H*H' ./ (attcount-1);
Sw=zeros(ecount);
for attidx=1:(attcount-1),
   Sw=Sw+squeeze(rH(:,attidx,:))*squeeze(rH(:,attidx,:))';
end
Sw=Sw./((attcount-1)*bootcount);

% code stolen from mikie to find fisher discriminant
D=Sb*pinv(Sw);
[eigVec,eigVal]=eig(D);
W=eigVec(:,1:(attcount-2));

eigmtx=ua(:,1:ecount);
Wp=eigmtx*W;
showkern(cat(3,Wp,zeros(size(Wp)),zeros(size(Wp))),kernfmt,iconside)

pH=H'*W;
peH=zeros(size(pH));
rH=squeeze(lH(1:ecount,:,2:end,:,:));
for attidx=1:(attcount-1),
   pH(attidx,:)=mean(squeeze(rH(:,attidx,:))'*W,1);
   peH(attidx,:)=std(squeeze(rH(:,attidx,:))'*W,1,1) .* ...
       sqrt((resampcount-1)/resampcount) .* sqrt(resampcount);
end

eigpatches=ua' * mpatches(:,1:attcount);
pP=eigpatches(1:ecount,2:end)'*W;
eigdistr=ua' * fpatches;
pD=eigdistr(1:ecount,:)'*W;

for ii=1:attcount-2,
   subplot(3,attcount-2,attcount-2+ii);
   errorbar(pH(:,ii),peH(:,ii));
   title(sprintf('dim %d kernel proj',ii));
   a=axis;
   axis([0 attcount a(3) a(4)]);
   
   subplot(3,attcount-2,2*(attcount-2)+ii);
   plot(0.5,pD(:,ii),'kx');
   hold on
   plot(pP(:,ii));
   hold off
   title(sprintf('dim %d target proj',ii));
   a=axis;
   axis([0 attcount a(3) a(4)]);
end




figure(3);
clf

subplot(4,1,1);
errorbar(1:attcount-1,c1mean(2:end),c1std(2:end));
hold on
plot([1 attcount-1],[c1mean(1) c1mean(1)],'k:');
hold off
title('DC response vs attidx');
ylabel('mean rate (Hz)');

subplot(4,1,2);
errorbar(1:attcount-1,a1mean(2:end),a1std(2:end));
hold on
plot([1 attcount-1],[a1mean(1) a1mean(1)],'k:');
hold off
title('global gain vs attidx');
ylabel('mod (x baseline)');

disp('not plotting prediction results because they aren''t done yet.');


return



hs=subplot(4,1,3);
paircount=(attcount-1)*(attcount-2);
pstr={};
pcount=0;
for p1=2:attcount-1,
   for p2=(p1+1):attcount,
      pcount=pcount+1;
      pstr{pcount}=sprintf('%d-%d',p1-1,p2-1);
   end
end
Tp=cat(4,pxc(1,1:end,:),pxcp(1,1:end,:,:));
Tp=reshape(Tp,size(Tp,2)*size(Tp,3),size(Tp,4))';
Txc=cat(4,xc(1,1:end,:),xcp(1,1:end,:,:));
Txc=reshape(Txc,size(Txc,2)*size(Txc,3),size(Txc,4))';

imagesc(Txc);
hold on
for ii=1:size(Tp,1),
   if ii==1,
      tt=find(Tp(ii,:)<PATTA);
   else
      tt=find(Tp(ii,:)<PATT);
   end
   plot(tt,ii.*ones(size(tt)),'kx');
end
hold off

axis xy
colorbar
title(sprintf('prediction correlations (pair p<%.3f all p<%.3f)',...
              PATT,PATTA));
set(hs,'YTick',1:length(pstr)+1,'YTickLabel',{'all',pstr{:}});
ylabel('pair');
set(hs,'XTickLabel',{'no','dc','dc/g','dc/l','no--','--','gl','loc'});
xlabel('model');

hs=subplot(4,1,4);

Tp=cat(4,pxc(1,2:end,:),pxcp(1,2:end,:,:));
Tp=reshape(Tp,size(Tp,2)*size(Tp,3),size(Tp,4))';

imagesc(-Tp,[-1 0]);
hold on
for ii=1:size(Tp,1),
   if ii==1,
      tt=find(Tp(ii,:)<PATTA);
   else
      tt=find(Tp(ii,:)<PATT);
   end
   plot(tt,ii.*ones(size(tt)),'kx');
end
hold off

%imagesc(Tp<PATTA,[0 1]);
axis xy
colorbar
title(sprintf('prediction significance (pair p<%.3f all p<%.3f)',...
              PATT,PATTA));
set(hs,'YTick',1:length(pstr)+1,'YTickLabel',{'all',pstr{:}});
ylabel('pair');
set(hs,'XTickLabel',{'dc','dc/g','dc/l','--','gl','loc'});
xlabel('model');

set(gcf,'PaperOrientation','Portrait','PaperPosition',[0.25 0.25 6 10.5]);
colormap(hot)

figure(3);
clf
for ii=1:length(targlist),
   subplot(length(targlist),1,ii);
   imagesc(bigpatches(:,:,targlist(ii)),[0 255]);
   title(sprintf('target %d (id %d)',ii,targlist(ii)));
   axis off
   axis image
end
colormap(gray);
set(gcf,'PaperOrientation','Portrait','PaperPosition',[0.25 0.25 8 10.5]);

if 0,
   %print the figures
   cd /auto/k1/david/docs/sfn02/CELLS/
   fprintf('printing to %s!\n',pwd);
   print('-depsc','-f1',sprintf('%s.kern.eps',cellid));
   print('-depsc','-f2',sprintf('%s.sig.eps',cellid));
   print('-depsc','-f3',sprintf('%s.targ.eps',cellid));
end
if 0,
   %print the figures
   fprintf('printing to printer\n');
   print -f1 -Pgcolor
   %print -f2
end



return

keyboard


figure(3);
clf
ps={'k-','g--','c--','r--','b--'};
ps2={'kx','ro','bo','go','co'};
labels={'A','1','2','3','4'};
ll=zeros(attcount,1);

spaceshow=1:20;

subplot(4,1,1);
attshowset=1:attcount;
Hs0=squeeze(std(Hpc(:,respidx,1,1,:),0,5)*sqrt(bootcount));
Hs0(find(Hs0==0))=1;
for attidx=attshowset,
   
   %a=errorbar(Hm(spaceshow,:),Hs(spaceshow,:),ps{attidx});
   a=plot(H(spaceshow,attidx)./Hs0(spaceshow),ps{attidx});
   %a=plot(H(spaceshow,attidx),ps{attidx});
   ll(attidx)=a(1);
   hold on
end
aa=axis;
axis([0 length(spaceshow)+1 aa(3) aa(4)]);
hold off
title(sprintf('%s PC domain STRF',cellid));
xlabel('pc idx');
legend(ll(find(ll)),labels{find(ll)});

hs=subplot(4,1,2);
paircount=(attcount-1)*(attcount-2);
Td=[];
pstr={};
pcount=0;
for p1=2:attcount-1,
   for p2=(p1+1):attcount,
      pcount=pcount+1;
      Td=[Td (abs(Hm(spaceshow,p1)-Hm(spaceshow,p2)) ./ ...
              sqrt(mean(Hs(spaceshow,[p1 p2]).^2,2)))];
      pstr{pcount}=sprintf('%d-%d',p1-1,p2-1);
   end
end

%imagesc(spaceuses,1:attdocount,pstrfatt(:,:,respidx),[0 1])
imagesc(Td',[0 3]);
axis xy
colorbar
set(hs,'YTickLabel',pstr);
ylabel('pair');
xlabel('PC');



keyboard

hold on
for pidx=1:attdocount,
   lidx=find(pstrfatt(pidx,:,respidx)<PATT);
   if length(lidx)>0,
      plot(lidx,pidx,'kx');
   end
end
hold off
ylabel('1-2  1-3  1-4  2-3  2-4  3-4');
xlabel(sprintf('signif pc idx (p<%.3f)',PATT));
set(hs,'YTickLabel',[]);
ttick=floor(str2num(get(hs,'XTickLabel')));
ttick(find(ttick==0))=1;
set(hs,'XTickLabel',spaceuses(ttick));
colormap(flipud(hot));
colorbar
drawnow

hs=subplot(3,1,3);
imagesc(attcc,[-1 1]);
title('strf/target cc');
xlabel('target');
ylabel('strf');
axis image

if ~isfield(z,'rbinned'),
   return
end

% display different kernels
figure(2);
clf
respidx=min([respcount 7]);
smark={'k','g','c','b','r'};
%rm=reshape(brmean(:,:,2:end),respcount,attcount-1,noisesets);
%rmall=squeeze(mean(rm(respidx,:,1)));
rmall=0;
rmin=min(min(min(min(rbinned(:,respidx,2:attcount,:)))))-rmall;
rmax=max(max(max(max(rbinned(:,respidx,2:attcount,:)))))-rmall;
spaceshow=spacelimp;
rowcount=ceil(spaceshow/5);
for eigidx=1:spaceshow,
   subplot(rowcount,ceil(spaceshow/rowcount),eigidx);
   
   if 0,
      imagesc(squeeze(rbinned(:,respidx,1:attcount,eigidx))',[rmin,rmax]);
      hold on
      for attidx=1:attcount,
         plot(mpbinned(attidx,eigidx),attidx,'kx');
         plot(mpbinned(attidx,eigidx),attidx,'wo');
      end
      hold off
      
      axis off;
   else
      hold on
      labels={'A','1','2','3','4'};
      ll=zeros(attcount,1);
      attshowset=[1 2 3 4 5];
         for attidx=attshowset,
            r=rbinned(:,respidx,attidx,eigidx);
            %r=rbinned(:,respidx,attidx,eigidx)-...
            %  mean(rbinned(:,respidx,attidx,eigidx));
            a=errorbar(sbin(:,eigidx),r,...
                       squeeze(rbinnedstd(:,respidx,attidx,eigidx)).*2,...
                       [smark{attidx},'-']);
            %set(a,'LineWidth',2);
            ll(attidx)=a(1);
            a=plot(sbin(mpbinned(attidx,eigidx),eigidx),...
                 r(mpbinned(attidx,eigidx)),'kx');
            set(a,'LineWidth',2);
         end
         hold off
         axis([sbin(1,eigidx),sbin(end,eigidx),rmin,rmax]);
         axis([sbin(1,eigidx),sbin(end,eigidx),rmin,rmax]);
         if eigidx==spaceshow,
            legend(ll(find(ll)),labels{find(ll)});
         end
   end
   title(sprintf('%s: pc%d r=%d',cellid,spaceusep(eigidx),respidx));
end
colormap(hot);
set(gcf,'PaperOrientation','Landscape','PaperPosition',[0.25 0.25 10.5 8]);
%set(gcf,'PaperOrientation','Landscape','PaperPosition',[0.25 2 10.5 4]);

figure(4);
clf

hs=subplot(3,1,1);
imagesc(perfatt0(:,:,respidx),[0 1]);
axis xy
hold on
for pidx=1:attdocount,
   lidx=find(perfatt0(pidx,:,respidx)<PATT);
   if length(lidx)>0,
      plot(lidx,pidx,'kx');
   end
end
hold off
ylabel('1-2  1-3  1-4  2-3  2-4  3-4');
xlabel(sprintf('signif pc idx (p<%.3f)',PATT));
title('erf0 diff');
set(hs,'YTickLabel',[]);
ttick=floor(str2num(get(hs,'XTickLabel')));
ttick(find(ttick==0))=1;
set(hs,'XTickLabel',spaceusep(ttick));
colormap(flipud(hot));
colorbar

hs=subplot(3,1,2);
imagesc(perfatt(:,:,respidx),[0 1]);
axis xy
hold on
for pidx=1:attdocount,
   lidx=find(perfatt(pidx,:,respidx)<PATT);
   if length(lidx)>0,
      plot(lidx,pidx,'kx');
   end
end
hold off
ylabel('1-2  1-3  1-4  2-3  2-4  3-4');
xlabel(sprintf('signif pc idx (p<%.3f)',PATT));
title('erf diff');
set(hs,'YTickLabel',[]);
ttick=floor(str2num(get(hs,'XTickLabel')));
ttick(find(ttick==0))=1;
set(hs,'XTickLabel',spaceusep(ttick));
colormap(flipud(hot));
colorbar

hs=subplot(3,1,3);
imagesc(perfatt1(:,:,respidx),[0 1]);
axis xy
hold on
for pidx=1:attdocount,
   lidx=find(perfatt1(pidx,:,respidx)<PATT);
   if length(lidx)>0,
      plot(lidx,pidx,'kx');
   end
end
hold off
ylabel('1-2  1-3  1-4  2-3  2-4  3-4');
xlabel(sprintf('signif pc idx (p<%.3f)',PATT));
title('erf diff individual PCs');
set(hs,'YTickLabel',[]);
ttick=floor(str2num(get(hs,'XTickLabel')));
ttick(find(ttick==0))=1;
set(hs,'XTickLabel',spaceusep(ttick));
colormap(flipud(hot));
colorbar
set(gcf,'PaperOrientation','Portrait','PaperPosition',[0.25 0.25 8 10.5]);


drawnow

return

figure(4);
clf
for eigidx=1:spacelimp,
   subplot(ceil(spacelimp/6),6,eigidx);
   pdata=[squeeze(perfatt0(:,eigidx,:)); ...
          squeeze(perfatt(:,eigidx,:))]';
   imagesc(pdata,[0 1]);
   
   hold on
   for respidx=1:respcount,
      tt=find(pdata(respidx,:)<PATT);
      if length(tt)>0,
            plot(tt,respidx,'bx');
      end
   end
   hold off
   
   if eigidx==1,
      title(sprintf('%s eigidx=%d',cellid,eigidx));
   else
      title(sprintf('eigidx=%d',eigidx));
   end
end
xlabel('1-2 1-3 1-4 2-3 2-4 3-4');
colormap(flipud(gray))
set(gcf,'PaperOrientation','Landscape','PaperPosition',[0.25 0.25 10.5 8]);


return

% switch to BW plotting scheme
torder=get(0,'DefaultAxesColorOrder');
set(0,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-|--|:|-.');
figure(4);
clf
subplot(3,1,1);
% find tuned dims
if min([pstrf(:); perf(:)])>PFILT,
   PFILT=min(pstrf(:));
end
pdata=[sum(sum(pstrf<PFILT,2),3) sum(sum(perf<PFILT,1),3)'];
plot(pdata);
legend('strf','erf');
title(sprintf('%s stimulus tuning',cellid));
ylabel(sprintf('n(att/resp) p<%.3f',PFILT));
xlabel('PC idx');

subplot(3,1,2);
% find att mod dims
pdata=[sum(sum(pstrfatt<PATT,1),3)' ...
       sum(sum(perfatt0<PATT,1),3)' ...
       sum(sum(perfatt<PATT,1),3)'];
plot(spaceusep,pdata);
legend('strf','erf0','erf');
title(sprintf('%s attention modulation',cellid));
ylabel(sprintf('n(att/resp) mod p<%.3f',PATT));
xlabel('PC idx');

subplot(3,1,3);
% find att mod dims
pdata=[sum(sum(pstrfatt(:,:,respidx)<PATT,1),3)' ...
       sum(sum(perfatt0(:,:,respidx)<PATT,1),3)' ...
       sum(sum(perfatt(:,:,respidx)<PATT,1),3)'];
plot(spaceusep,pdata);
legend('strf','erf0','erf');
title(sprintf('%s attention modulation respidx=%d',cellid,respidx));
ylabel(sprintf('n(att/resp) mod p<%.3f',PATT));
xlabel('PC idx');

% return to default (color) plotting scheme
set(0,'DefaultAxesColorOrder',torder);


% code to dump all cells

sql='SELECT * FROM sRunData WHERE batch=11 ORDER BY cellid';
rundata=mysql(sql);
for ii=1:length(rundata),
   kerncomp4res(rundata(ii).cellid,11);
end






