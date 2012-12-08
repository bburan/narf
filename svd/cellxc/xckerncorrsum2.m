
batches=[23; 29; 32];
batches2=[86; 83; 85; 84; 87; 88];

dbopen;

res=loadbigdata(batches);
res2=loadbigdata(batches2);
predxc=res.predxc;
prederr=res.prederr;
predxcst=res2.predxc;
prederrst=res2.prederr;
predp=res.predp;
celllist=res.celllist;

batchcount=length(batches);
cellcount=length(celllist);

batchcount=length(batches);
cellcount=length(celllist);

badcellidx=find(strcmp(celllist,'r0150B') | strcmp(celllist,'93G83A'));
predxc(badcellidx,1,:,:,:)=nan;
predxc(badcellidx,:,:,1,:)=nan;

expstr={'Rev','GR','NR'};
expidx1=[1 1 2 ];
expidx2=[2 3 3 ];
cnfidx1=[1 1 2 ];
cnfidx2=[2 2 3 ];
nloutidx1=[2 2 2 ];

stcnfidx1=[1 1 1 ];
stcnfidx2=[1 1 1 ];
posnlidx1=[1 1 2 ];
posnlidx2=[2 3 3 ];
stnlidx1=[4 4 2 ];
stnlidx2=[2 3 3 ];

compcount=length(expidx1);

%
% summarize prediction stuff
%
figure(1);
clf

for ii=1:compcount,
   classexists=~isnan(predxc(:,expidx1(ii))) & ~isnan(predxc(:,expidx2(ii)));
   allclassexists=sum(~isnan(predxc(:,:,1)),2)==3;
   
   p1=predxc(:,expidx1(ii),1,cnfidx1(ii),nloutidx1(ii));
   p2=predxc(:,expidx2(ii),1,cnfidx2(ii),nloutidx1(ii));
   e1=prederr(:,expidx1(ii),1,cnfidx1(ii),nloutidx1(ii));
   e2=prederr(:,expidx2(ii),1,cnfidx2(ii),nloutidx1(ii));
   
   %withingoodidx=find(classexists & p1-e1>0 & p2-e2>0);
   withingoodidx=find(classexists);
   allwithingoodidx=find(classexists & p1-e1>0 & p2-e2>0);
   %allwithingoodidx=find(allclassexists & p1-e1>0 & p2-e2>0);
   %allwithingoodidx=find(allclassexists);
   
   pr1=predxc(:,expidx1(ii),1,1,nloutidx1(ii));
   pr2=predxc(:,expidx2(ii),1,1,nloutidx1(ii));
   pst1=predxcst(:,1,1,stcnfidx1(ii),stnlidx1(ii));
   pst2=predxcst(:,1,1,stcnfidx2(ii),stnlidx2(ii));
   pstgr=predxcst(:,3+expidx2(ii),1,stcnfidx2(ii),4);
   ppos1=predxcst(:,2,1,stcnfidx1(ii),posnlidx1(ii));
   ppos2=predxcst(:,2,1,stcnfidx2(ii),posnlidx2(ii));
   %pneg1=predxcst(:,3,1,stcnfidx1(ii),posnlidx1(ii));
   %pneg2=predxcst(:,3,1,stcnfidx2(ii),posnlidx2(ii));
   pneg1=zeros(size(ppos1));
   pneg2=zeros(size(ppos1));
   
   pall=[pr1 pst1 ppos1 pneg1 ppos2 pst2 pstgr pr2];
   
   subplot(compcount,2,(ii-1)*2+1);
   %m=median(pall(withingoodidx,:));
   m=nanmean(pall(withingoodidx,:));
   bar(m);
   title(sprintf('%s v %s (n=%d)',expstr{expidx1(ii)},...
                 expstr{expidx2(ii)},length(withingoodidx)));
   axis([0 length(m)+1 0 0.5]);
   
   subplot(compcount,2,(ii-1)*2+2);
   %m=median(pall(allwithingoodidx,:));
   m=nanmean(pall(allwithingoodidx,:));
   bar(m);
   title(sprintf('%s v %s w/in good (n=%d)',expstr{expidx1(ii)},...
                 expstr{expidx2(ii)},length(allwithingoodidx)));
   axis([0 length(m)+1 0 0.5]);   
end
colormap(gray);
fullpage('portrait');

figure(2);
clf

classexists=~isnan(predxc(:,1)) & ~isnan(predxc(:,2));
withingoodidx=find(classexists);

subplot(2,2,1);

set1=predxcst(withingoodidx,5,1,1,4);
set2=predxcst(withingoodidx,1,1,1,4);
e1=prederrst(withingoodidx,5,1,1,4);
e2=prederrst(withingoodidx,1,1,1,4);

set1hi=find(set1-e1 > set2+e2);
set2hi=find(set1+e1 < set2-e2);
sigidx=ones(size(set1))*2;
sigidx(set1hi)=3;
sigidx(set2hi)=1;

n1='gratrev';
n2='review';

[p,m]=plotcomp(set1,set2,n1,n2,[-0.4 1.0 -0.4 1.0],sigidx);

subplot(2,2,2);

set1=predxcst(withingoodidx,1,1,1,2);
set2=predxcst(withingoodidx,1,1,1,4);
e1=prederrst(withingoodidx,1,1,1,2);
e2=prederrst(withingoodidx,1,1,1,4);

set1hi=find(set1-e1 > set2+e2);
set2hi=find(set1+e1 < set2-e2);
sigidx=ones(size(set1))*2;
sigidx(set1hi)=3;
sigidx(set2hi)=1;

n1='gratrev rt';
n2='review rt';

[p,m]=plotcomp(set1,set2,n1,n2,[-0.4 1.0 -0.4 1.0],sigidx);

subplot(2,2,3);

set1=predxcst(withingoodidx,2,1,1,2);
set2=predxcst(withingoodidx,2,1,1,1);
e1=prederrst(withingoodidx,2,1,1,2);
e2=prederrst(withingoodidx,2,1,1,1);

set1hi=find(set1-e1 > set2+e2);
set2hi=find(set1+e1 < set2-e2);
sigidx=ones(size(set1))*2;
sigidx(set1hi)=3;
sigidx(set2hi)=1;

n1='gratrev pos';
n2='review pos';

[p,m]=plotcomp(set1,set2,n1,n2,[-0.4 1.0 -0.4 1.0],sigidx);

subplot(2,2,4);
pr=cat(3,predxc(withingoodidx,1:2,1,1:2,2),...
       predxcst(withingoodidx,[1 5],1,1:2,4));
pr=permute(pr,[1 4 2 3]);
mpr=reshape(mean(pr),4,2);

bar(mpr);
axis square
title('within-/across class sep,non-sep');

colormap(gray);

drawnow

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

clear tunecomp spaceresp

goodidx2=1:cellcount;

batchtrans=[4 2 3];

for ii=1:cellcount,
   % load kern data for corresponding cell and fill
   sql=['SELECT * FROM sRunData WHERE cellid="',celllist{ii},'"',...
        ' AND batch=',num2str(batches2(1))];
   rundata=mysql(sql);
   
   if ~ismember(ii,badcellidx) & sum(~isnan(predxc(ii,:,1,1,1)))>1 ...
         & ismember(ii,goodidx2) & length(rundata)>=1,
      
      resfile=[rundata(1).respath,rundata(1).resfile,'.gz'];
      fprintf('%d: %s\n',ii,resfile);
      r0=zload(resfile);
      
      bm=find(r0.goodbatch(batchtrans));
      
      sidx=zeros(size(bm));
      for jj=1:length(bm),
         sidx(jj)=find(r0.goodbatchrange==batchtrans(bm(jj)));
         
         % tweak--force gratrev kernels to look unbiased
         if bm(jj)==3,
            r0.strf(sidx(jj)).powunbiased(:)=1;
         end
      end
      
      strf=r0.strf(sidx);
      
      r=compkernels(strf);
      
      timecorr(bm,bm,ii,1)=r.tc;
      etimecorr(bm,bm,ii,1)=r.etc;
      spacecorr(bm,bm,ii,1)=r.sc;
      espacecorr(bm,bm,ii,1)=r.esc;
      spacecorr(bm,bm,ii,2)=r.scp;
      spacecorr(bm,bm,ii,3)=r.scn;
      
      spacecorr2(bm,bm,ii,:)=reshape(r0.estxc(sidx,sidx,:),...
                                     length(bm),length(bm),1,4);
      
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

expstr={'Rev','GR','NR'};
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





