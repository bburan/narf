% function z=kencompres(cellid,batch)
%
% load results from kernfile in sRunData and display attentional
% modulation info
%
% z=data from kernfile
%
function z=kerncompres(runidx,batch)

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
attdocount=z.attdocount;
spacelims=z.spacelims;
spacelimp=z.spacelimp;
respcount=z.respcount;
PFILT=z.PFILT;
PATT=z.PATT;

respidx=min([z.respcount,7]);

% draw the strf one eigenvector at a time, and their sum
ecount=min([10 z.spacelims]);
spacecount=z.spacecount;
eigmtx=z.ua(:,z.spaceuses(1:ecount));
eigtargets=eigmtx'*(z.mpatches-repmat(z.bsmean',[1 5]));

eigmtx=repmat(eigmtx,[1 1 attcount]);
targmtx=reshape(z.mpatches,spacecount,1,attcount);
if strcmp(z.kernfmt,'pfft') | strcmp(z.kernfmt,'pfftgr'),
   targmtx(37,:)=0;
end
for attidx=1:attcount,
   for jj=1:ecount,
      eigmtx(:,jj,attidx)=eigmtx(:,jj,attidx).*...
          z.Hpc(jj,respidx,attidx,1);
   end
   %targmtx(:,:,attidx)=targmtx(:,:,attidx)-targmtx(:,:,1);
end

for attidx=1:attcount,
   targmtx(:,:,attidx)=targmtx(:,:,attidx) ./ ...
       max(abs(targmtx(:,:,attidx))) .* max(max(max(abs(eigmtx(:,:,2:end)))));
end
strf=sum(eigmtx,2);
eigmtx=cat(2,targmtx,eigmtx,strf);
%for attidx=1:attcount,
%  eigmtx(:,:,attidx)=eigmtx(:,:,attidx)./max(max(abs(eigmtx(:,:,attidx))));
%end

attcc=zeros(attcount);
for xx=1:attcount,
   for yy=1:attcount,
      attcc(xx,yy)=xcov(strf(:,xx),targmtx(:,yy),0,'coeff');
   end
end


figure(3);
clf
kernfmt=z.kernfmt;
if strcmp(kernfmt,'pfftgr'),
   kernfmt='pfft';
   z.iconside=[sqrt(spacecount*2) sqrt(spacecount*2)];
end
showkern(eigmtx,kernfmt,z.iconside,{},0);
for eigidx=0:ecount,
   subplot(attcount,ecount+2,eigidx+1);
   if eigidx==0,
      title(sprintf('%s r=%d',cellid,respidx));
   else
      title(sprintf('pc%d',z.spaceuses(eigidx)));
   end
   xlabel('');
end
subplot(attcount,ecount+2,ecount+2);
title('STRF');

for attidx=1:attcount,
   subplot(attcount,ecount+2,attidx*(ecount+2)-ecount);
   ylabel(sprintf('attidx=%d',attidx));
end

figure(1);
clf
ps={'k-','r--','b--','g--','c--'};
ps2={'kx','ro','bo','go','co'};
labels={'A','1','2','3','4'};
ll=zeros(attcount,1);

subplot(3,1,1);
hold on
attshowset=1:attcount;
for attidx=attshowset,
   a=errorbar(z.Hpcfull(z.spaceuses,respidx,attidx,1),...
              std(z.Hpcfull(:,respidx,attidx,2:end),0,4), ...
              ps{attidx});
   ll(attidx)=a(1);
   plot(z.Hpc(:,respidx,attidx,1),ps2{attidx});
end
aa=axis;
axis([0 size(z.Hpcfull,1)+1 aa(3) aa(4)]);
hold off
title(sprintf('%s PC domain STRF',cellid));
xlabel('pc idx');
legend(ll(find(ll)),labels{find(ll)});

hs=subplot(3,1,2);
%imagesc(z.spaceuses,1:attdocount,z.pstrfatt(:,:,respidx),[0 1])
imagesc(z.pstrfatt(:,:,respidx),[0 1]);
axis xy
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


