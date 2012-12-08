% function kerncompsum(batchid,latidx)
%
function kerncompsum(batchid,latidx)

dbopen;

if ~exist('batchid'),
   batchid=11;
end
if ~exist('latidx'),
   latidx=1;
end

sql=['SELECT sResults.*',...
     ' FROM sResults INNER JOIN sRunData',...
     ' ON sResults.runid=sRunData.id',...
     ' WHERE sResults.batch=',num2str(batchid),...
     ' ORDER BY cellid'];
resdata=mysql(sql);

for ii=1:length(resdata),
   eval(char(resdata(ii).matstr));
end

okidx=zeros(length(preddata),1);
tt=0;
for ii=1:length(preddata),
   if ~isempty(preddata(ii).pxc),
      okidx(ii)=1;
      tt=tt+1;
      kerndata(tt)=preddata(ii);
   end
end
ttt=find(okidx)';

cellids={kerndata.cellid};
[sc,scidx]=sort(cellids);
scidx=1:length(cellids);

if 0,
   keyboard
end

PFILT=0.01;
PATT=0.008; % 0.008;
PXCT=0.05;
fprintf('batch: %d  latidx=%d  tuning: pstrf/perf\n',batchid,latidx);
xccount=zeros(length(kerndata),3,2);
predxc=zeros(length(kerndata),3,4,2);
predxccross=zeros(length(kerndata),3,4,5,2);
nondampedcount=zeros(length(kerndata),5);
ncount=zeros(length(kerndata),5);
plist=ones(length(kerndata),3,2);
pxctick=0;
fprintf('CELL: (in/out) noatt dcatt globalatt fullatt p: pdc pglob pfull\n');
for ii=1:length(kerndata),
   kidx=scidx(ii);
   fprintf('%-5s:',kerndata(kidx).cellid);
   
   if isfield(kerndata(kidx),'pxc'),
      predcount=size(kerndata(kidx).predxc,2);
      for predidx=1:predcount,
         fprintf('%5.2f/%5.2f ',kerndata(kidx).predxc(1,predidx),...
                 kerndata(kidx).predxc(2,predidx));
      end
      
      xcc=size(kerndata(kidx).pxc,3);
      fprintf('p<%.2f %.2f %.2f',kerndata(kidx).pxc(1,2:4,1));
      fprintf('  %.2f %.2f %.2f',kerndata(kidx).pxc(1,2:4,2));
      plist(ii,:,1:xcc)=kerndata(kidx).pxc(1,2:4,1:xcc);
      
      xcc=size(kerndata(kidx).predxc,3);
      predxc(ii,:,:,1:xcc)=...
          reshape(kerndata(kidx).predxc(:,:,1:xcc),1,3,4,xcc);
      xccount(ii,:,:)=(plist(ii,:,:)<PXCT);
      acc=size(kerndata(kidx).predxccross,3);
      predxccross(ii,:,:,1:acc,1:xcc)=...
          reshape(kerndata(kidx).predxccross(:,:,1:acc,1:xcc),1,3,4,acc,xcc);
      fprintf(' %d%d%d %d%d%d',xccount(ii,:));
      pxctick=pxctick+1;
      ndl=length(kerndata(kidx).nondampedcount);
      if ndl>0 & ~strcmp(kerndata(kidx).cellid,'m0000'),
         nondampedcount(kidx,1:ndl)=kerndata(kidx).nondampedcount(:)';
      end
      ndl=length(kerndata(kidx).ncount);
      if ndl>0 & ~strcmp(kerndata(kidx).cellid,'m0000'),
         ncount(kidx,1:ndl)=kerndata(kidx).ncount(:)';
      end
   else
      attcount=size(kerndata(kidx).pstrf,3);
      spacelim=size(kerndata(kidx).perf,2);
      for attidx=1:attcount,
         perfcount=sum(kerndata(kidx).perf(attidx,:,latidx)<PFILT);
         pstrfcount=sum(kerndata(kidx).pstrf(:,1,attidx)<PFILT);
         fprintf('%-2d/%-2d ',pstrfcount,perfcount);
      end
      
   end
   
   fprintf('\n');
end

% skip m0000 on summary numbers:
realidx=1:(size(predxc,1)-1);

% only look at local attention for cells with good snr
localgoodidx=find(sum(nondampedcount(:,2:end)>0,2)>=1);

ii=realidx';
nid=zeros(length(ii),2,3,2);
STT=2;

for xcidx=1:2,
   % test for cells with signif DC mod
   
   % stderr on randomized preds
   st=predxc(ii,3,2,1);
   st(find(st==0))=1;
   
   % significant shifts only if predxc significantly above randomized
   % AND significantly above zero
   nid(:,1,1,xcidx)=(plist(ii,1,xcidx)<PXCT & abs(predxc(ii,1,2,1))./st>STT);
   nid(:,2,1,xcidx)=(abs(predxc(ii,1,2,1))./st>STT);
   
   % test for cells with signif GLOBAL mod
   st=predxc(ii,3,3,1);
   st(find(st==0))=1;
   nid(:,1,2,xcidx)=(plist(ii,2,xcidx)<PXCT & abs(predxc(ii,1,3,1))./st>STT);
   nid(:,2,2,xcidx)=(abs(predxc(ii,1,3,1))./st>STT);
   
   % test for cells with signif LOCAL mod
   st=predxc(ii,3,4,1);
   st(find(st==0))=1;
   nid(:,1,3,xcidx)=(plist(ii,3,xcidx)<PXCT & abs(predxc(ii,1,4,1))./st>STT);
   nid(:,2,3,xcidx)=(abs(predxc(ii,1,4,1))./st>STT);
end

nsum=squeeze(sum(nid));

fprintf('MODEL (predsig/N tunesig/Ngood) mean/rand meandiff mediandiff\n');
fprintf('DC (%d/%d %d/%d) %5.3f/%5.3f %5.3f %5.3f\n',nsum(:,1,:),...
        mean(predxc(realidx,1:2,2,1)),mean(predxc(realidx,1,2,1)-predxc(realidx,2,2,1)),...
        median(predxc(realidx,1,2,1)-predxc(realidx,2,2,1)));
fprintf('GL (%d/%d %d/%d) %5.3f/%5.3f %5.3f %5.3f\n',nsum(:,2,:),...
        mean(predxc(realidx,1:2,3,1)),mean(predxc(realidx,1,3,1)-predxc(realidx,2,3,1)),...
        median(predxc(realidx,1,3,1)-predxc(realidx,2,3,1)));
fprintf('LO (%d/%d %d/%d) %5.3f/%5.3f %5.3f %5.3f\n',nsum(:,3,:),...
        mean(predxc(realidx,1:2,4,1)),mean(predxc(realidx,1,4,1)-predxc(realidx,2,4,1)),...
        median(predxc(realidx,1,4,1)-predxc(realidx,2,4,1)));


%disp('skipping plot');
%return

figure(1)
clf
catidx=zeros(max(realidx),1);
arange=[-0.15 0.15 0 30];

subplot(3,3,1);
catidx(realidx)=1;
catidx(find(xccount(realidx,1)))=2;
%plotcomp(predxc(realidx,2,2,1),predxc(realidx,1,2,1),...
%         'no att','dc',[0 1 0 1],catidx,{'k','r'},{'filled','filled'});
aa=find(~nid(:,1,1) & nid(:,2,1));
bb=find(nid(:,1,1));
histcomp(predxc(aa,1,2,1)-predxc(aa,2,2,1),...
         predxc(bb,1,2,1)-predxc(bb,2,2,1),...
         'no att','dc','xc',arange);
fprintf('dc good mean diff= %.3f\n',...
        mean(predxc(bb,1,2,1)-predxc(bb,2,2,1)));

subplot(3,3,2);
catidx(realidx)=1;
catidx(find(xccount(realidx,2)))=2;
%plotcomp(predxc(realidx,2,3,1),predxc(realidx,1,3,2),...
%         'no att','global',[0 1 0 1],catidx,{'k','r'},{'filled','filled'});
aa=find(~nid(:,1,2) & nid(:,2,2));
bb=find(nid(:,1,2));
histcomp(predxc(aa,1,3,1)-predxc(aa,2,3,1),...
         predxc(bb,1,3,1)-predxc(bb,2,3,1),...
         'no att','global','xc',arange);
if length(bb)>0,
   fprintf('global good mean diff= %.3f\n',...
           mean(predxc(bb,1,3,1)-predxc(bb,2,3,1)));
else
   fprintf('no global diff\n');
end

subplot(3,3,3);
aa=find(~nid(:,1,3) & nid(:,2,3));
bb=find(nid(:,1,3));
histcomp(predxc(aa,1,4,1)-predxc(aa,2,4,1),...
         predxc(bb,1,4,1)-predxc(bb,2,4,1),...
         'no att','global','xc',arange);
if length(bb)>0,
   fprintf('local good mean diff= %.3f\n',...
           mean(predxc(bb,1,4,1)-predxc(bb,2,4,1)));
else
   fprintf('no local diff\n');
end


if 0, 
subplot(3,2,3);
%plotcomp(predxc(realidx,1,3,2)-predxc(realidx,2,1,2),...
%         predxc(realidx,1,4,2)-predxc(realidx,2,1,2),...
%         'global-dc','local-dc',[-0.5 0.5 -0.5 0.5],catidx);
aa=find(~nid(:,1,1) & ~nid(:,1,2) & ~nid(:,1,3) & nid(:,2,1));
bb=find(nid(:,1,1) | nid(:,1,2) | nid(:,1,3));
histcomp(predxc(aa,1,1,1),predxc(bb,1,1,1),...
         'no att','','xc',[0 1 arange(3:4)]);

subplot(3,3,5);
%catidx(realidx)=1;
%catidx(find(xccount(realidx,5)))=2;
%plotcomp(predxc(realidx,2,3,2),predxc(realidx,1,3,2),...
%         'dc att','global',[0 1 0 1],catidx,{'k','r'},{'filled','filled'});
aa=find(~nid(:,1,2) & nid(:,2,2));
bb=find(nid(:,1,2));
histcomp(predxc(aa,1,3,2)-predxc(aa,2,3,2),...
         predxc(bb,1,3,2)-predxc(bb,2,3,2),...
         'dc att','global','xc',arange);

subplot(3,3,6);
%catidx(realidx)=1;
%catidx(find(xccount(realidx,6)))=2;
%plotcomp(predxc(realidx,2,4,2),predxc(realidx,1,4,2),...
%         'dc att','local',[0 1 0 1],catidx,{'k','r'},{'filled','filled'});
aa=find(~nid(:,1,3) & nid(:,2,3));
bb=find(nid(:,1,3));
histcomp(predxc(aa,1,4,2)-predxc(aa,2,4,2),...
         predxc(bb,1,4,2)-predxc(bb,2,4,2),...
         'dc att','local','xc',arange);
end

hs=subplot(3,3,4);
counts=[sum(plist(ii,3,1)<PXCT & nid(:,2,3)) ...
        sum(plist(ii,2,1)<PXCT & nid(:,2,2) & ...
            (plist(ii,3,1)>=PXCT | ~nid(:,2,3))) ...
        sum(plist(ii,1,1)<PXCT & nid(:,2,1) & ...
            (plist(ii,2,1)>=PXCT | ~nid(:,2,2)) & ...
            (plist(ii,3,1)>=PXCT | ~nid(:,2,3))) ...
        sum(plist(ii,1,1)>=PXCT & nid(:,2,1) & ...
            (plist(ii,2,1)>=PXCT | ~nid(:,2,2)) & ...
            (plist(ii,3,1)>=PXCT | ~nid(:,2,3)))];
pie(counts+0.001,...
    {['local ',num2str(counts(1))],['global ',num2str(counts(2))],...
     ['DC ',num2str(counts(3))],['None ',num2str(counts(4))]});
title('significant pred summ');

hs=subplot(3,3,5);
ii=realidx';
counts=cat(3,[sum(plist(ii,:,1)<PXCT,2)==0 & nid(:,2,1)  ...
            zeros(length(ii),2)],...
           [plist(ii,1,1)<PXCT & nid(:,2,1) ...
            zeros(length(ii),2)],...
           [plist(ii,2,1)<PXCT & plist(ii,1,1)<PXCT & nid(:,2,2)...
            plist(ii,2,1)<PXCT & plist(ii,1,1)>=PXCT & nid(:,2,2) ...
            zeros(length(ii),1)],...
           [plist(ii,3,1)<PXCT & plist(ii,1,1)<PXCT & nid(:,2,3) ...
            plist(ii,3,1)<PXCT & plist(ii,1,1)>=PXCT & nid(:,2,3) ...
            zeros(length(ii),1)]);
counts=squeeze(sum(counts,1))';
bar(counts,'stacked');
set(hs , 'XTickLabel',{'None','DC','global','local'});
title('signif pred summ 2');

hs=subplot(3,3,6);
counts=[sum(plist(ii,3,2)<PXCT & nid(:,2,3)) ...
        sum(plist(ii,2,2)<PXCT & nid(:,2,2) & ...
            (plist(ii,3,2)>=PXCT | ~nid(:,2,3))) ...
        sum(plist(ii,1,2)<PXCT & nid(:,2,1) & ...
            (plist(ii,2,2)>=PXCT | ~nid(:,2,2)) & ...
            (plist(ii,3,2)>=PXCT | ~nid(:,2,3))) ...
        sum(plist(ii,1,2)>=PXCT & nid(:,2,1) & ...
            (plist(ii,2,2)>=PXCT | ~nid(:,2,2)) & ...
            (plist(ii,3,2)>=PXCT | ~nid(:,2,3)))];
pie(counts+0.001,...
    {['local ',num2str(counts(1))],['global ',num2str(counts(2))],...
     ['DC ',num2str(counts(3))],['None ',num2str(counts(4))]});
title('significant tuning summ');

hs=subplot(3,3,7);
%plotcomp(plist(ii,1,1),plist(ii,1,2),'DC pred p','tune p',...
%         [0 1 0 1],2-nid(:,1,1,1));
hist(plist(ii,1,1),40);
title('DC pred p');

hs=subplot(3,3,8);
%plotcomp(plist(ii,2,1),plist(ii,2,2),'GL pred p','tune p',...
%         [0 1 0 1],2-nid(:,1,2,1));
hist(plist(ii,2,1),40);
title('GL pred p');

hs=subplot(3,3,9);
%plotcomp(plist(ii,3,1),plist(ii,3,2),'LO pred p','tune p',...
%         [0 1 0 1],2-nid(:,1,3,1));
hist(plist(ii,3,1),40);
title('LO pred p');


set(gcf,'PaperOrientation','portrait','PaperPosition',[0.25 0.25 8 10.5]);

return




[squeeze(predxc(ii,:,4,2)) plist(ii,3,2) ...
 abs(predxc(ii,2,4,2))./predxc(ii,3,4,2)]

figure(2)
hs=subplot(3,1,1);
hist(plist(:,:,1))

hs=subplot(3,1,2);
hist(plist(:,:,2))


vv=[1 0; -0.5 sqrt(3)./2; -0.5 -sqrt(3)/2];
xx=squeeze(predxc(localgoodidx,1,2:4,1)) * vv;


fprintf('batch: %d  latidx=%d  attmod: pstrf(S)/perf(E)/perf0(R)\n',...
        batchid,latidx);
bmod=zeros(length(kerndata),3);
strfbar=zeros(length(kerndata),20);
erfbar=zeros(length(kerndata),20);
s1=' SER';
for ii=1:length(kerndata),
   kidx=scidx(ii);
   attcount=size(kerndata(kidx).pstrfatt,1);
   spacelim=size(kerndata(kidx).perfatt,2);
   fprintf('%-5s (%-2d): ',kerndata(kidx).cellid,spacelim);
   for attidx=1:attcount,
      ee=min([size(kerndata(kidx).pstrfatt,2)]);
      pstrfcount=sum(kerndata(kidx).pstrfatt(attidx,ee,1)<PATT);
      if size(kerndata(kidx).perfatt,2)>0,
         ee=min([size(kerndata(kidx).perfatt,2) 10]);
         %perfcount=sum(kerndata(kidx).perfatt(attidx,end,latidx)<PATT);
         perfcount=sum(kerndata(kidx).perfatt(attidx,ee,latidx)<PATT);
      else
         perfcount=0;
      end
      perf0count=sum(kerndata(kidx).perfatt0(attidx,end,latidx)<PATT);
      fprintf('%-2d/%-2d/%-2d ',pstrfcount,perfcount,perf0count);
      if pstrfcount>0,
         bmod(kidx,1)=1;
      end
      if perfcount>0,
         bmod(kidx,2)=1;
      end
      if perf0count>0,
         bmod(kidx,3)=1;
      end
   end
   
   % thing to see where you get max - spacelim=26???
   strflim=size(kerndata(kidx).pstrfatt,2);
   strfbar(ii,1:strflim)=sum((kerndata(kidx).pstrfatt(:,:,1)<PATT),1)>0;
   strfbar(ii,strflim+1:end)=strfbar(ii,strflim);
   erflim=size(kerndata(kidx).perfatt,2);
   erfbar(ii,1:erflim)=sum((kerndata(kidx).perfatt(:,:,1)<PATT),1)>0;
   erfbar(ii,erflim+1:end)=erfbar(ii,erflim);
   
   fprintf('%s%s%s\n',s1(bmod(kidx,1)+1),s1(bmod(kidx,2)*2+1),...
           s1(bmod(kidx,3)*3+1));
end
fprintf('Mod counts: %d/%d/%d\n',sum(bmod,1));

fprintf(' ecount:');
fprintf(' %2d',1:size(strfbar,2));
fprintf('\n');
fprintf('strfbar:');
fprintf(' %2d',sum(strfbar,1));
fprintf('\n erfbar:');
fprintf(' %2d',sum(erfbar,1));
fprintf('\n');

fprintf('pxc count: %d %d %d/%d\n',sum(xccount,1),pxctick);

return


cellid={preddata(ttt).cellid};
cellcount=length(cellid);
cellcount=cellcount+1;
cellid{cellcount}='mean';

runidx=find(okidx);
runidx(cellcount)=0;

predxc=cat(5,preddata(ttt).predxc);
tpredxc=predxc;
tpredxc(find(isnan(tpredxc)))=0;
predxc=cat(5,predxc,mean(tpredxc,5));
%attsimxc=cat(4,preddata(ttt).attsimxc);
%attsimxc=cat(4,attsimxc,mean(attsimxc,4));

TSON=0;
if TSON,
   targsim=cat(4,preddata(ttt).targsim);
   targsim=cat(4,targsim,mean(targsim,4));
end

% make the last entry in each summary statistic the mean across
% cells
if length(preddata(ttt(1)).p)>1,
   p=repmat(zeros(size(preddata(ttt(1)).p)),[1 1 1 length(ttt)]);
   for ii=1:length(ttt),
      pcount0=min([size(preddata(ttt(ii)).p,1) size(p,1)]);
      pcount=min([size(preddata(ttt(ii)).p,2) size(p,2)]);
      pcount2=min([size(preddata(ttt(ii)).p,3) size(p,3)]);
      p(1:pcount0,1:pcount,1:pcount2,ii)=...
          preddata(ttt(ii)).p(1:pcount0,1:pcount,1:pcount2);
   end
   p=cat(4,p,mean(p,4));
else
   p=zeros(size(predxc));
end

% display prediction scores for this latency bin (4 means 0-100 ms after
% fixation onset)
latcount=size(predxc,3);
latbase=min([latcount 1])-1;
if ~exist('latidx'),
   latidx=min([latcount 6]);
end
PCLOW=10;
if ~exist('pccount'),
   pccount=18
end
nlidx=1;

% PTHRESH=0.012; % p < 1-(1-0.05)^(1/4)
PTHRESH=0.008; % p < 1-(1-0.05)^(1/6
%PTHRESH=0.05
figure(1);
clf
minx=min(min(min(predxc(:,:,latidx,:))));
maxx=max(max(max(predxc(:,:,latidx,:))));
if minx<-maxx & maxx>0,
   minx=-maxx;
end
for ii=1:cellcount,
   subplot(ceil(cellcount/4),4,ii);
   imagesc(predxc(:,:,latidx,nlidx,ii),[minx maxx]);
   axis image
   axis off
   title(sprintf('%s (%d)',cellid{ii},runidx(ii)));
end
colorbar
colormap(hot);

attcount=size(predxc,1);
attpcount=size(p,1);

disp('cell  SMT (pred in/out) (att comp 1-2 1-3 1-4 2-3 2-4 3-4)');
pxc=zeros(cellcount,attcount);
pxcout=zeros(cellcount,attcount);
pxcatt=[];
sigcount=[0 0];
attpairs=[1 2; 1 3; 1 4; 2 3; 2 4; 3 4];

xcallin=[];
xcallout=[];
pall=[];

attbetter=zeros(cellcount-1,attcount-1);
for ii=1:cellcount,
   fprintf('%-5s',cellid{ii});
   for attidx=1:attcount,
      %tpredxc=predxc(:,:,9,nlidx,ii);
      tpredxc=predxc(:,:,latidx,nlidx,ii);
      %tpredxc(find(isnan(tpredxc)))=0;
      pxc(ii,attidx)=tpredxc(attidx,attidx);
      pxcout(ii,attidx)=mean(tpredxc(attidx,[2:attidx-1 (attidx+1):end]));
      
      if attidx>1 & ii<cellcount,
         xcallin=[xcallin tpredxc(attidx,attidx).*ones(1,attcount-2)];
         xcallout=[xcallout tpredxc(attidx,[2:attidx-1 (attidx+1):end])];
      end
      
      pxcout(ii,attidx)=mean(tpredxc(attidx,[2:attidx-1 (attidx+1):end]));
      fprintf(' %5.2f/%5.2f',pxc(ii,attidx),pxcout(ii,attidx));
      
      if ii<cellcount & attidx>1 & pxc(ii,attidx)>pxcout(ii,attidx),
         attbetter(ii,attidx-1)=1;
      end
   end
   
   % print out in vs out pred xcs. transfer over to pairwise
   % comparison if it looks like it's gonna help.
   
   if preddata(ttt(1)).p~=0,
      for attidx=1:attpcount,
         pall=[pall; p(attidx,pccount,latidx-latbase,ii)];
         if p(attidx,pccount,latidx-latbase,ii)<=PTHRESH,
            sextra='*';
            %pxcatt[pxcatt; pxc(ii,attidx) pxcout(ii,attidx)];
            pxcatt=[pxcatt; pxc(ii,attpairs(attidx,1)+1) ...
                    pxc(ii,attpairs(attidx,2)+1)];
         else
            sextra=' ';
         end
         fprintf('%s%5.3f',sextra,p(attidx,pccount,latidx-latbase,ii));
      end
      
      % if qualifies, mark as cell with significant attentional effect
      if min(p(1:attpcount,PCLOW,latidx-latbase,ii))<=PTHRESH,
         fprintf('*%d',PCLOW);
         sigcount(2)=sigcount(2)+1;
      end
      if min(p(1:attpcount,pccount,latidx-latbase,ii))<=PTHRESH,
         fprintf('*%d',pccount);
         sigcount(1)=sigcount(1)+1;
      end
   else
      pxcatt=[0 0];
   end
   fprintf('\n');
end
fprintf(['sig att: %d/%d %.2f%% (%dPC) %d/%d %.2f%% (%dPC)',...
         ' p<%.3f respidx=%d\n'],...
        sigcount(2),cellcount-1,sigcount(2)./(cellcount-1).*100,...
        PCLOW,...
        sigcount(1),cellcount-1,sigcount(1)./(cellcount-1).*100,...
        pccount,PTHRESH,latidx);
xcallin=xcallin(find(~isnan(xcallout)));
xcallout=xcallout(find(~isnan(xcallout)));
[pp,mm,ss]=randpairtest(xcallin,xcallout,2000);

fprintf('mean xc in/out: %.3f/%.3f (p<%.3f): %d/%d (%.1f%%) attstates\n',...
        mean(xcallin),mean(xcallout),pp,...
        sum(xcallin>xcallout),length(xcallout),...
        sum(xcallin>xcallout)/length(xcallout)*100);

figure(2);
clf
subplot(2,3,1);
scatter(xcallout(:),xcallin(:));
hold on
plot([0 1.0],[0 1.0],'k--');
plot(nanmean(xcallout),nanmean(xcallin),'r+','LineWidth',2);
hold off
axis equal
xlabel('pcx out');
ylabel('pcx in');
title('attention preds in vs. out');

subplot(2,3,2);
scatter(mean(pxcout(:,2:end),2),pxc(:,1));
hold on
plot([0 1.0],[0 1.0],'k--');
plot(median(mean(pxcout(:,2:end),2)),median(pxc(:,1)),...
     'r+','LineWidth',2);
hold off
axis equal
xlabel('pcx att in');
ylabel('pcx all att');
title('data lim: all vs. atten preds');

subplot(2,3,3);
mxcin=zeros(latcount,1);
mxcout=zeros(latcount,1);
for attidx=2:attcount,
   mxcin=mxcin+squeeze(predxc(attidx,attidx,:,nlidx,end)) ./ (attcount-1);
   mxcout=mxcout+squeeze(sum(predxc(attidx,[2:attidx-1 (attidx+1):end],...
                                     :,nlidx,end),2)) ./ ...
          ((attcount-2)*(attcount-1));
end
plot([mxcout mxcin]);
legend('out','in');
%scatter(pxcatt(:,2),pxcatt(:,1));
%hold on
%plot([0 1.0],[0 1.0],'k--');
%plot(median(pxcatt(:,2)),median(pxcatt(:,1)),...
%     'r+','LineWidth',2);
%hold off
%axis equal
%xlabel('xc att out');
%ylabel('xc att in');
title(sprintf('att mod preds (latidx=%d,pccount=%d)',latidx,pccount));

subplot(2,3,4);
%plot(squeeze(sum((min(p(2:end,:,:,:),[],1)<=PTHRESH),4)),'LineWidth',2);
plot(squeeze(sum((min(p(1:end,:,:,:),[],1)<=PTHRESH),4)),'LineWidth',2);
hold on
%plot(squeeze(sum((min(p(2:end,:,latidx-latbase,:),[],1)<=PTHRESH),4)),...
%     'kx','LineWidth',2);
plot(squeeze(sum((min(p(1:end,:,latidx-latbase,:),[],1)<=PTHRESH),4)),...
     'kx','LineWidth',2);
hold off
xlabel('pc count');
ylabel('N signif atten');
legend('lat1','lat2');
title('sig att cells vs. pc count');

subplot(2,3,5);
hist(pall,15);
title('pair p values')

subplot(2,3,6);
hist(xcallin(:)-xcallout(:),15);
title('pair pred diff (in - out)');

if TSON,
   subplot(2,3,5);
   ts1=targsim(:,pccount,1,:);
   ts2=targsim(:,pccount,2,:);
   tp=p(:,pccount,latidx-latbase,:);
   sigsim=[ts1(:) ts2(:) tp(:)];
   scatter(sigsim(:,1),sigsim(:,3));
   title('att diff vs targ diff (unnorm)');
   xlabel('targ1.*targ2');
   ylabel('att p');
   
   subplot(2,3,6);
   scatter(sigsim(find(tp<PTHRESH),2),sigsim(find(tp<PTHRESH),3));
   title('att diff vs targ diff (norm)');
   xlabel('angle(targ1,targ2)');
   ylabel('att p');
end

%keyboard
return

figure(3);
clf
atxc=squeeze(attsimxc(:,:,pccount,:));
minx=min(atxc(:));
maxx=max(atxc(:));
if minx<-maxx & maxx>0,
   minx=-maxx;
end
for ii=1:cellcount-1,
   subplot(ceil(cellcount/4),4,ii);
   imagesc(atxc(:,:,ii),[minx maxx]);
   axis image
   axis off
   title(sprintf('atxc %s (%d)',cellid{ii},runidx(ii)));
end
colorbar
colormap(hot);



disp('cell  T   (pred in/out)');
for ii=1:cellcount,
   fprintf('%-5s',cellid{ii});
   for attidx=1:attcount,
      pxc=predxc(attidx,attidx,2,ii);
      pxcout=mean(predxc(attidx,[2:attidx-1 (attidx+1):end],latidx,ii));
      fprintf(' %5.3f (%4.1f/%4.1f)',p(attidx,1,1,ii),pxc,pxcout);
   end
   
   % if qualifies, mark as cell with significant attentional effect
   if min(p(:,1,1,ii))<PTHRESH,
      fprintf('*');
   end
   fprintf('\n');
end


