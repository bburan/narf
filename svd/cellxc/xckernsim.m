
batches=[78; 81];
nlidx=[3; 3];
cnfidx=[3; 2];
expstr={'NR','GR'};

dbopen;
sql=['SELECT DISTINCT cellid FROM sRunData WHERE batch=',...
     num2str(batches(1)),' AND cellid<"r0306" ORDER BY cellid'];
rundata=mysql(sql);

celllist={rundata.cellid};
badcells=(strcmp(celllist,'r0150B') | strcmp(celllist,'93G83A') | ...
          strcmp(celllist,'modelhic') | strcmp(celllist,'modelhis') ...
          );

cellcount=length(rundata);
batchcount=length(batches);

%
% pull out features from kernels
%
timecorr=nan.*zeros(batchcount,batchcount,cellcount,3);
etimecorr=nan.*zeros(batchcount,batchcount,cellcount,3);
spacecorr=nan.*zeros(batchcount,batchcount,cellcount,3);
espacecorr=nan.*zeros(batchcount,batchcount,cellcount,3);
spacecorr2=nan.*zeros(batchcount,batchcount,cellcount,4);
adaptidx=nan.*zeros(cellcount,batchcount);
eadaptidx=nan.*zeros(cellcount,batchcount);
suppidx=nan.*zeros(cellcount,batchcount);
esuppidx=nan.*zeros(cellcount,batchcount);
sepidx=nan.*zeros(cellcount,batchcount);
esepidx=nan.*zeros(cellcount,batchcount);
tempresp=nan.*zeros(14,batchcount,cellcount);
predxc=nan.*zeros(cellcount,batchcount);
prederr=nan.*zeros(cellcount,batchcount);

clear tunecomp spaceresp
for ii=1:cellcount,
   clear strf
   
   for batchidx=1:batchcount,
      % load kern data for corresponding cell and fill
      sql=['SELECT * FROM sRunData WHERE cellid="',celllist{ii},'"',...
           ' AND batch=',num2str(batches(batchidx))];
      rundata=mysql(sql);
      
      if length(rundata)>0,
         resfile=[rundata.respath,rundata.resfile,'.gz'];
      else
         resfile=[];
      end
      if isempty(resfile) | ~exist(resfile,'file'),
         fprintf('cell %s batch %d absent\n',celllist{ii},...
                 batches(batchidx));
         res(ii).bad=1;
      else    
         r0=zload(resfile);
         
         strf(batchidx,1)=r0.strf(nlidx(batchidx));
         strf(batchidx,1).parms.sfscount=r0.params.sfscount;
         strf(batchidx,1).parms.sfsstep=r0.params.sfsstep;
         
         predxc(ii,batchidx)=r0.predxc(cnfidx(batchidx),nlidx(batchidx));
         prederr(ii,batchidx)=r0.prederr(cnfidx(batchidx),nlidx(batchidx));
      end
   end
   
   if sum(isnan(predxc(ii,:)))==0, % ie, both exist
      r=compkernels(strf);
      bm=1:2;
      timecorr(bm,bm,ii,1)=r.tc;
      etimecorr(bm,bm,ii,1)=r.etc;
      spacecorr(bm,bm,ii,1)=r.sc;
      espacecorr(bm,bm,ii,1)=r.esc;
      spacecorr(bm,bm,ii,2)=r.scp;
      spacecorr(bm,bm,ii,3)=r.scn;
      spacecorr(bm,bm,ii,4)=r.sc0;
      spacecorr(bm,bm,ii,5)=r.scp0;
      spacecorr(bm,bm,ii,6)=r.scn0;
      
      adaptidx(ii,bm)=r.adaptidx;
      eadaptidx(ii,bm)=r.eadaptidx;
      suppidx(ii,bm)=r.suppidx;
      esuppidx(ii,bm)=r.esuppidx;
      sepidx(ii,bm)=r.sepidx;
      esepidx(ii,bm)=r.esepidx;
      
      tempresp(:,bm,ii)=r.tempresp;
      
      if ~exist('spaceresp','var'),
         spaceresp=nan.*zeros(size(r.hspace,1),batchcount,cellcount);
      end
      spaceresp(1:size(r.hspace,1),bm,ii)=r.hspace;
   end
end

figure(1);
clf

for ii=1:2,
   subplot(3,3,ii);
   
   set1=squeeze(tempresp(:,ii,:));
   goodidx=find(~isnan(predxc(:,1)+predxc(:,2)));
   set1=set1(:,goodidx);
  
   % normalize temporal responses to have equal weight
   for jj=1:size(set1,2),
      set1(:,jj)=set1(:,jj)./norm(set1(:,jj));
   end
   
   set2=cumsum(set1,1);
   
   plot([7 (size(set1,1)*14-7)],[0 0],'k--');
   hold on
   plot(7:14:(size(set1,1)*14-7),mean(set1,2),'k-','LineWidth',2);
   %plot(7:14:(size(set1,1)*14-7),mean(set2,2),'--');
   hold off
   
   title(sprintf('%s mean tempresp (n=%d)',...
                 expstr{ii},size(set1,2)));
   axis square
   a=axis;
   axis([a(1) a(2) -0.25 0.75]);
   
   xlabel('time lag (ms)');
end


subplot(3,3,3);
set1=adaptidx(goodidx,1);
set2=adaptidx(goodidx,2);

p1=predxc(goodidx,1);
p2=predxc(goodidx,2);
e1=prederr(goodidx,1);
e2=prederr(goodidx,2);

catidx=ones(length(goodidx),1).*1;
%catidx(find(p1-e1<0 | p2-e2<0))=2;

plotcomp(set1,set2,[expstr{1},' adapt'],expstr{2},...
         [0 1 0 1],catidx);

set1=squeeze(spacecorr(1,2,goodidx,:));
for jj=1:6,
   
   subplot(3,3,3+jj);
   
   d1=set1(:,jj);
   d2=[];
   
   sunits='xc';
   axisrange=[-0.5 1 0 15];
   
   [p1,p2]=histcomp(d1(:),d2(:),'a','b',...
                       sunits,axisrange);
end


keyboard


expstr={'NR','GR'};
cnfstr=expstr;
outstr={'lin','thr','ethr','full'};

expidx1=[1 1 3 ];
expidx2=[2 3 2 ];
stidx1=[4 4 2 ];
stidx2=[2 3 3 ];
cnfidx1=[1 1 3 ];
cnfidx2=[2 3 2 ];
nloutidx1=[2 2 2 ];

DOSPACE=1;

rowcount=length(expidx1);
colcount=3;

SEC=2.0;

fprintf('DOSPACE=%d SEC=%.1f\n',DOSPACE,SEC);

% figure out some stuff about how to group cells
revandgridx=find(~isnan(predxc(:,1)) & ~isnan(predxc(:,2)));
goodidx4=find(~isnan(predxc(:,1)) & ~isnan(predxc(:,2)) & ...
              ~isnan(predxc(:,3)));
%goodidx2=goodidx3;

colcount=5;

figure(4);
clf
for ii=1:rowcount,
   
   subplot(rowcount,colcount,(ii-1)*colcount+1);
   
   set1=squeeze(tempresp(:,ii,:));
   goodidx=find(~isnan(set1(1,:)));
   goodidx=intersect(goodidx,revandgridx);
   set1=set1(:,goodidx);
   
   m1=mean(set1([1 2 end-1 end],:));
   set1=set1-repmat(m1,size(set1,1),1);
   
   v1=max(abs(set1),[],1);
   set1=set1./repmat(v1,size(set1,1),1);
   set2=cumsum(set1,1);
   
   plot([7 (size(set1,1)*14-7)],[0 0],'k--');
   hold on
   plot(7:14:(size(set1,1)*14-7),mean(set1,2));
   plot(7:14:(size(set1,1)*14-7),mean(set2,2),'--');
   hold off
   
   title(sprintf('%s mean tempresp (n=%d)',...
                 expstr{ii},size(set1,2)));
   axis square
   a=axis;
   axis([a(1) a(2) -0.5 1.5]);
   
   if ii==rowcount,
      xlabel('time lag (ms)');
   end
   
   subplot(rowcount,colcount,(ii-1)*colcount+2);
   set1=adaptidx(:,expidx1(ii));
   set2=adaptidx(:,expidx2(ii));
   goodidx=find(~isnan(set1) & ~isnan(set2));
   goodidx=intersect(goodidx,goodidx2);
   set1=set1(goodidx);
   set2=set2(goodidx);
   
   p1=predxc(goodidx,expidx1(ii),1,cnfidx1(ii),nloutidx1(ii));
   p2=predxc(goodidx,expidx2(ii),1,cnfidx2(ii),nloutidx1(ii));
   e1=prederr(goodidx,expidx1(ii),1,cnfidx1(ii),nloutidx1(ii));
   e2=prederr(goodidx,expidx2(ii),1,cnfidx2(ii),nloutidx1(ii));
   
   catidx=ones(length(goodidx),1).*1;
   catidx(find(p1-e1<0 | p2-e2<0))=2;
   
   plotcomp(set1,set2,[expstr{expidx1(ii)},' adapt'],expstr{expidx2(ii)},...
            [0 1 0 1],catidx);
   
   % get correlation between various componets of the spatial (or
   % temporal) response function
   if DOSPACE==0,
      % pull out temp corr data
      set1=squeeze(timecorr(expidx1(ii),expidx2(ii),:,:));
      set2=squeeze(timecorr(expidx1(ii),expidx2(ii),:,:)+...
                   etimecorr(expidx1(ii),expidx2(ii),:,:).*SEC);
      nroot=sprintf('%s/%s time',...
                    expstr{expidx1(ii)},expstr{expidx2(ii)});
   elseif DOSPACE==1,
      % pull out spatial corr data
      set1=squeeze(spacecorr(expidx1(ii),expidx2(ii),:,:));
      set2=squeeze(spacecorr(expidx1(ii),expidx2(ii),:,:)+...
                   espacecorr(expidx1(ii),expidx2(ii),:,:).*SEC);
      nroot=sprintf('%s/%s space',...
                    expstr{expidx1(ii)},expstr{expidx2(ii)});
   elseif DOSPACE==2,
      % pull out spatial corr data
      set1=squeeze(spacecorr2(expidx1(ii),expidx2(ii),:,1:3));
      set2=squeeze(spacecorr2(expidx1(ii),expidx2(ii),:,1:3));
      nroot=sprintf('%s/%s space2',...
                    expstr{expidx1(ii)},expstr{expidx2(ii)});
   end
   nastr={[nroot,' sig'],[nroot,' pos sig'],[nroot,' neg sig']};
   nbstr={'ns','ns','ns'};
   
   set1(find(set1>0.999))=0.999;
   cumspacecomp(:,:,ii)=set1(goodidx4,:);
   
   set1=set1(goodidx,:);
   set2=set2(goodidx,:);
   set1(find(isnan(set1)))=0;
   
   % classify as either sig or ns based on w/in class pred accuracy
   p1=predxc(goodidx,expidx1(ii),1,cnfidx1(ii),nloutidx1(ii));
   p2=predxc(goodidx,expidx2(ii),1,cnfidx2(ii),nloutidx1(ii));
   e1=prederr(goodidx,expidx1(ii),1,cnfidx1(ii),nloutidx1(ii));
   e2=prederr(goodidx,expidx2(ii),1,cnfidx2(ii),nloutidx1(ii));
   
   % classify as to review st-sep kernel is better than hybrid at
   % rev pred
   pst1=predxcst(goodidx,1,1,1,stidx1(ii));
   pst2=predxcst(goodidx,1,1,1,stidx2(ii));
   est1=prederrst(goodidx,1,1,1,stidx1(ii));
   est2=prederrst(goodidx,1,1,1,stidx2(ii));
   
   cat1idx=find(p1-e1>=0 & p2-e2>=0 & pst1>=pst2);
   cat2idx=find(p1-e1<0 | p2-e2<0);
   cat3idx=find(p1-e1>=0 & p2-e2>=0 & pst1<pst2);
   %cat1idx=find(pst1-est1-est2>pst2);
   %cat2idx=find(pst1-est1-est2<=pst2);
   %cat3idx=[];
   
   fprintf('1: %s vs 2: %s\n',...
           expstr{expidx1(ii)},expstr{expidx2(ii)});
   fprintf('cat1: both good within, pred1>pred2 (n=%d)\n',length(cat1idx));
   for cc=1:length(cat1idx),
      tidx=cat1idx(cc);
      fprintf('%6s: %6.3f v %6.3f : %6.3f  %6.3f  %6.3f\n',...
              celllist{goodidx(tidx)},pst1(tidx),pst2(tidx),...
                 set1(tidx,:));
   end
   fprintf('cat2: at least one bad w/in (n=%d)\n',length(cat2idx));
   for cc=1:length(cat2idx),
      tidx=cat2idx(cc);
      fprintf('%6s: %6.3f v %6.3f : %6.3f  %6.3f  %6.3f\n',...
              celllist{goodidx(tidx)},pst1(tidx),pst2(tidx),...
                 set1(tidx,:));
   end
   fprintf('cat3: both good within, pred1<pred2 (n=%d)\n',length(cat3idx));
   for cc=1:length(cat3idx),
      tidx=cat3idx(cc);
      fprintf('%6s: %6.3f v %6.3f : %6.3f  %6.3f  %6.3f\n',...
              celllist{goodidx(tidx)},pst1(tidx),pst2(tidx),...
                 set1(tidx,:));
   end
   
   for jj=1:3,
      
      % classify as to whether spatial correlations are
      % significantly different
      %d1=set1(find(set2(:,jj)<1),jj);
      %d2=set1(find(set2(:,jj)>=1),jj);
      
      % divide up good and bad predictors
      d1=set1([cat1idx;cat3idx],jj);
      d2=set1([cat2idx],jj);
      
      % get rid of negative values?
      %d1=abs(d1);
      %d2=abs(d2);
      
      % don't divide up good and bad predictors
      %d1=set1([cat1idx;cat2idx;cat3idx],jj);
      %d2=[];
      
      % classify as all-three or not
      %goodidx3=find(sum(isnan(predxc(goodidx,:)),2)==0);
      %d1=set1(goodidx3,jj);
      %d2=set1(setdiff(1:size(set1,1),goodidx3),jj);
      
      subplot(rowcount,colcount,(ii-1)*colcount+2+jj);
      
      sunits='xc';
      axisrange=[-0.5 1 0 15];
      
      %seta=[p1-p2];
      %setb=[set1(:,jj)];
      %catidx=ones(size(setb));
      %catidx(find(p1-e1<0 | p2-e2<0 | pst1<=pst2))=2;
      %
      %plotcomp(seta,setb,[expstr{expidx1(ii)},nastr{jj}],...
      %         expstr{expidx2(ii)},[-1 1 -1 1],catidx);
      
      [p1,p2]=histcomp(d1(:),d2(:),nastr{jj},nbstr{jj},...
                       sunits,axisrange);
   end
end

colormap(gray);
fullpage('landscape')


figure(2);
clf

subplot(2,2,1);
batchcomp([92;92],[3;1],3,0,0);

subplot(2,2,2);
batchcomp([94;94],[3;1],3,0,0);

subplot(2,2,3);
batchcomp([92;92],[3;1],3,0,2);

subplot(2,2,4);
batchcomp([94;94],[3;1],3,0,2);
fullpage('portrait');

