% function z=kencomp4res(cellid,batch)
%
% load results from kernfile in sRunData and display attentional
% modulation info ... kernfile generated from kerncomp4
%
% z=data from kernfile
%
% created SVD 10/18/02 - hacked from kerncompres.m
%
function z=kerncomp4res(runidx,batch)

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

rcsetstrings;

if length(rundata)>1,
   disp('cellres-kerncomp: only displaying first batch');
end

fprintf('Loading %s...\n',[rundata(1).respath,rundata(1).kernfile,'.gz']);
z=zload([rundata(1).respath,rundata(1).kernfile,'.gz']);

attcount=z.attcount;
bootcount=z.bootcount;
spacelims=z.spacelims;
respcount=z.respcount;
spacecount=z.spacecount;
resampcount=z.params.resampcount;
PFILT=z.PFILT;
PATT=z.PATT;
PATTA=z.PATTA;

respidx=min([z.respcount,7]);

% draw the strf one eigenvector at a time, and their sum
ecount=max(z.predspace); % min([40 spacelims]);
fprintf('ecount=%d\n',ecount);
erange=z.spaceuses(1:ecount);
eigtargets=z.ua(:,z.spaceuses(1:ecount))' * ...
    (z.mpatches-repmat(z.bsmean',[1 size(z.mpatches,2)]));

if isfield(z,'DAMPFACTOR'),
   DAMPFACTOR=z.DAMPFACTOR;
else
   DAMPFACTOR=1.0;
end
noiseidx=1;

Hm=squeeze(mean(z.lH(:,respidx,:,noiseidx,:),5));
Hmed=squeeze(median(z.lH(:,respidx,:,noiseidx,:),5));
Hs=squeeze(std(z.lH(:,respidx,:,noiseidx,:),0,5).*sqrt(bootcount));
%Hdamp=(abs(Hm)./(Hs.*DAMPFACTOR)).^2;
%Hdamp(find(Hdamp>1))=1;
Hdamp=(abs(Hm)./(Hs.*DAMPFACTOR));
Hdamp=(1-Hdamp.^(-2));
Hdamp=Hdamp.*(Hdamp>0);
Hdamp(find(isnan(Hdamp)))=0;

if 0,
   % damp out according to shrinkage filter
   H=Hm.*repmat(Hdamp(:,1),[1 attcount]);
elseif z.params.dodamping
   % damp out according to shrinkage filter
   fprintf('sfsidx max=%d\n',z.ddmax);
   H=Hm.*z.Hdamp;
else
   H=Hm;
end

a1mean=mean(z.a1(respidx,:,1,:,1),4);
a1std=std(z.a1(respidx,:,1,:,1),0,4).*sqrt(bootcount);
%rateconv=z.respfilterparms{2}(respidx)-z.respfilterparms{1}(respidx)
% scale c1 by 1000 to convert to Hz from spikes/ms
c1mean=mean(z.c1(respidx,:,1,:,1),4)*1000;
c1std=std(z.c1(respidx,:,1,:,1),0,4).*sqrt(bootcount)*1000;

kernmtx=z.ua(:,erange)*H(erange,:);
globmtx=kernmtx(:,1)*a1mean;
targmtx=reshape(z.mpatches(:,1:attcount),spacecount,1,attcount);
if strcmp(z.params.kernfmt,'pfft') | strcmp(z.params.kernfmt,'pfftgr'),
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
kernfmt=z.params.kernfmt;
if strcmp(kernfmt,'pfftgr'),
   kernfmt='pfft';
   z.iconside=[sqrt(spacecount*2) sqrt(spacecount*2)];
end
showkern(eigmtx,kernfmt,z.iconside,{},0);

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
rH=squeeze(z.lH(1:ecount,:,2:end,1,:));
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

eigmtx=z.ua(:,1:ecount);
Wp=eigmtx*W;
showkern(cat(3,Wp,zeros(size(Wp)),zeros(size(Wp))),kernfmt,z.iconside)

pH=H'*W;
peH=zeros(size(pH));
rH=squeeze(z.lH(1:ecount,:,2:end,:,:));
for attidx=1:(attcount-1),
   pH(attidx,:)=mean(squeeze(rH(:,attidx,:))'*W,1);
   peH(attidx,:)=std(squeeze(rH(:,attidx,:))'*W,1,1) .* ...
       sqrt((resampcount-1)/resampcount) .* sqrt(resampcount);
end

eigpatches=z.ua' * z.mpatches(:,1:attcount);
pP=eigpatches(1:ecount,2:end)'*W;
eigdistr=z.ua' * z.fpatches;
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
Tp=cat(4,z.pxc(1,1:end,:),z.pxcp(1,1:end,:,:));
Tp=reshape(Tp,size(Tp,2)*size(Tp,3),size(Tp,4))';
Txc=cat(4,z.xc(1,1:end,:),z.xcp(1,1:end,:,:));
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

Tp=cat(4,z.pxc(1,2:end,:),z.pxcp(1,2:end,:,:));
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
for ii=1:length(z.targlist),
   subplot(length(z.targlist),1,ii);
   imagesc(z.bigpatches(:,:,z.targlist(ii)),[0 255]);
   title(sprintf('target %d (id %d)',ii,z.targlist(ii)));
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
Hs0=squeeze(std(z.Hpc(:,respidx,1,1,:),0,5)*sqrt(bootcount));
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

%imagesc(z.spaceuses,1:attdocount,z.pstrfatt(:,:,respidx),[0 1])
imagesc(Td',[0 3]);
axis xy
colorbar
set(hs,'YTickLabel',pstr);
ylabel('pair');
xlabel('PC');



keyboard

hold on
for pidx=1:attdocount,
   lidx=find(z.pstrfatt(pidx,:,respidx)<PATT);
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
set(hs,'XTickLabel',z.spaceuses(ttick));
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
respidx=min([z.respcount 7]);
smark={'k','g','c','b','r'};
%rm=reshape(z.brmean(:,:,2:end),respcount,attcount-1,noisesets);
%rmall=squeeze(mean(rm(respidx,:,1)));
rmall=0;
rmin=min(min(min(min(z.rbinned(:,respidx,2:attcount,:)))))-rmall;
rmax=max(max(max(max(z.rbinned(:,respidx,2:attcount,:)))))-rmall;
spaceshow=z.spacelimp;
rowcount=ceil(spaceshow/5);
for eigidx=1:spaceshow,
   subplot(rowcount,ceil(spaceshow/rowcount),eigidx);
   
   if 0,
      imagesc(squeeze(z.rbinned(:,respidx,1:attcount,eigidx))',[rmin,rmax]);
      hold on
      for attidx=1:z.attcount,
         plot(z.mpbinned(attidx,eigidx),attidx,'kx');
         plot(z.mpbinned(attidx,eigidx),attidx,'wo');
      end
      hold off
      
      axis off;
   else
      hold on
      labels={'A','1','2','3','4'};
      ll=zeros(attcount,1);
      attshowset=[1 2 3 4 5];
         for attidx=attshowset,
            r=z.rbinned(:,respidx,attidx,eigidx);
            %r=z.rbinned(:,respidx,attidx,eigidx)-...
            %  mean(z.rbinned(:,respidx,attidx,eigidx));
            a=errorbar(z.sbin(:,eigidx),r,...
                       squeeze(z.rbinnedstd(:,respidx,attidx,eigidx)).*2,...
                       [smark{attidx},'-']);
            %set(a,'LineWidth',2);
            ll(attidx)=a(1);
            a=plot(z.sbin(z.mpbinned(attidx,eigidx),eigidx),...
                 r(z.mpbinned(attidx,eigidx)),'kx');
            set(a,'LineWidth',2);
         end
         hold off
         axis([z.sbin(1,eigidx),z.sbin(end,eigidx),rmin,rmax]);
         axis([z.sbin(1,eigidx),z.sbin(end,eigidx),rmin,rmax]);
         if eigidx==spaceshow,
            legend(ll(find(ll)),labels{find(ll)});
         end
   end
   title(sprintf('%s: pc%d r=%d',cellid,z.spaceusep(eigidx),respidx));
end
colormap(hot);
set(gcf,'PaperOrientation','Landscape','PaperPosition',[0.25 0.25 10.5 8]);
%set(gcf,'PaperOrientation','Landscape','PaperPosition',[0.25 2 10.5 4]);

figure(4);
clf

hs=subplot(3,1,1);
imagesc(z.perfatt0(:,:,respidx),[0 1]);
axis xy
hold on
for pidx=1:attdocount,
   lidx=find(z.perfatt0(pidx,:,respidx)<PATT);
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
set(hs,'XTickLabel',z.spaceusep(ttick));
colormap(flipud(hot));
colorbar

hs=subplot(3,1,2);
imagesc(z.perfatt(:,:,respidx),[0 1]);
axis xy
hold on
for pidx=1:attdocount,
   lidx=find(z.perfatt(pidx,:,respidx)<PATT);
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
set(hs,'XTickLabel',z.spaceusep(ttick));
colormap(flipud(hot));
colorbar

hs=subplot(3,1,3);
imagesc(z.perfatt1(:,:,respidx),[0 1]);
axis xy
hold on
for pidx=1:attdocount,
   lidx=find(z.perfatt1(pidx,:,respidx)<PATT);
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
set(hs,'XTickLabel',z.spaceusep(ttick));
colormap(flipud(hot));
colorbar
set(gcf,'PaperOrientation','Portrait','PaperPosition',[0.25 0.25 8 10.5]);


drawnow

return

figure(4);
clf
for eigidx=1:spacelimp,
   subplot(ceil(spacelimp/6),6,eigidx);
   pdata=[squeeze(z.perfatt0(:,eigidx,:)); ...
          squeeze(z.perfatt(:,eigidx,:))]';
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
if min([z.pstrf(:); z.perf(:)])>PFILT,
   PFILT=min(z.pstrf(:));
end
pdata=[sum(sum(z.pstrf<PFILT,2),3) sum(sum(z.perf<PFILT,1),3)'];
plot(pdata);
legend('strf','erf');
title(sprintf('%s stimulus tuning',cellid));
ylabel(sprintf('n(att/resp) p<%.3f',PFILT));
xlabel('PC idx');

subplot(3,1,2);
% find att mod dims
pdata=[sum(sum(z.pstrfatt<PATT,1),3)' ...
       sum(sum(z.perfatt0<PATT,1),3)' ...
       sum(sum(z.perfatt<PATT,1),3)'];
plot(z.spaceusep,pdata);
legend('strf','erf0','erf');
title(sprintf('%s attention modulation',cellid));
ylabel(sprintf('n(att/resp) mod p<%.3f',PATT));
xlabel('PC idx');

subplot(3,1,3);
% find att mod dims
pdata=[sum(sum(z.pstrfatt(:,:,respidx)<PATT,1),3)' ...
       sum(sum(z.perfatt0(:,:,respidx)<PATT,1),3)' ...
       sum(sum(z.perfatt(:,:,respidx)<PATT,1),3)'];
plot(z.spaceusep,pdata);
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






