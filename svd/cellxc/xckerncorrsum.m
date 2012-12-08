
SFFT=1;
if SFFT,
   batches=[23; 32; 29];
   batches2=[52; 60; 56];
   batches3=[86; 88; 87];
else
   batches=[26; 31; 28];
   batches2=[51; 59; 55];
   batches3=[83; 85; 84];
end

dbopen;

res=loadbigdata(batches);
res3=loadbigdata(batches3);
predxc=res.predxc;
prederr=res.prederr;
predxcst=res3.predxc;
prederrst=res3.prederr;
predp=res.predp;
celllist=res.celllist;

batchcount=length(batches);
cellcount=length(celllist);

badcellidx=find(strcmp(celllist,'r0150B') | strcmp(celllist,'93G83A'));
predxc(badcellidx,1,:,:,:)=nan;

timecorr=nan.*zeros(batchcount,batchcount,cellcount,3);
etimecorr=nan.*zeros(batchcount,batchcount,cellcount,3);
spacecorr=nan.*zeros(batchcount,batchcount,cellcount,3);
espacecorr=nan.*zeros(batchcount,batchcount,cellcount,3);
spacecorr2=nan.*zeros(batchcount,batchcount,cellcount,3);
espacecorr2=nan.*zeros(batchcount,batchcount,cellcount,3);
spacecorr3=nan.*zeros(batchcount,batchcount,cellcount,3);
espacecorr3=nan.*zeros(batchcount,batchcount,cellcount,3);
tunecorr=nan.*zeros(batchcount,batchcount,cellcount,3);
etunecorr=nan.*zeros(batchcount,batchcount,cellcount,3);
adaptidx=nan.*zeros(cellcount,batchcount);
eadaptidx=nan.*zeros(cellcount,batchcount);
sepidx=nan.*zeros(cellcount,batchcount);
esepidx=nan.*zeros(cellcount,batchcount);
%suppidx=nan.*zeros(cellcount,batchcount,2);
%esuppidx=nan.*zeros(cellcount,batchcount,2);

suppidx=nan.*zeros(14,batchcount,cellcount).*nan;
esuppidx=nan.*zeros(14,batchcount,cellcount).*nan;
tempresp=nan.*zeros(14,batchcount,cellcount);

clear tunecomp;

%goodidx2=find(~isnan(predxc(:,1)) & ~isnan(predxc(:,2)) & ...
%             ~isnan(predxc(:,3)));
goodidx2=1:cellcount;

for ii=1:cellcount,
   % load kern data for corresponding cell and fill
   if ~ismember(ii,badcellidx) & sum(~isnan(predxc(ii,:,1,1,1)))>1 ...
         & ismember(ii,goodidx2),
      
      r=xccomptune(celllist{ii},batches2);

      %close
      
      %tunecomp(ii)=r;
      
      bm=zeros(size(r.batch));
      for jj=1:length(r.batch),
         bm(jj)=find(r.batch(jj)==batches2);
      end
      
      timecorr(bm,bm,ii,1)=r.tsc;
      etimecorr(bm,bm,ii,1)=r.etsc;
      timecorr(bm,bm,ii,2)=r.tpc;
      etimecorr(bm,bm,ii,2)=r.etpc;
      timecorr(bm,bm,ii,3)=r.tnc;
      etimecorr(bm,bm,ii,3)=r.etnc;
      
      spacecorr(bm,bm,ii,1)=r.ssc;
      espacecorr(bm,bm,ii,1)=r.essc;
      spacecorr(bm,bm,ii,2)=r.ssposc;
      espacecorr(bm,bm,ii,2)=r.essposc;
      spacecorr(bm,bm,ii,3)=r.ssnegc;
      espacecorr(bm,bm,ii,3)=r.essnegc;
      
      % sum over all timelags
      spacecorr2(bm,bm,ii,1)=r.sac;
      espacecorr2(bm,bm,ii,1)=r.esac;
      spacecorr2(bm,bm,ii,2)=r.spc;
      espacecorr2(bm,bm,ii,2)=r.espc;
      spacecorr2(bm,bm,ii,3)=r.snc;
      espacecorr2(bm,bm,ii,3)=r.esnc;
      
      % sum up to peak time lag
      spacecorr3(bm,bm,ii,1)=r.spc;
      espacecorr3(bm,bm,ii,1)=r.espc;
      spacecorr3(bm,bm,ii,2)=r.spposc;
      espacecorr3(bm,bm,ii,2)=r.espposc;
      spacecorr3(bm,bm,ii,3)=r.spnegc;
      espacecorr3(bm,bm,ii,3)=r.espnegc;
      
      adaptidx(ii,bm)=r.adaptidx;
      eadaptidx(ii,bm)=r.eadaptidx;
      sepidx(ii,bm)=r.sepidx;
      esepidx(ii,bm)=r.esepidx;
      %suppidx(ii,bm,1)=r.suppidx;
      %esuppidx(ii,bm,1)=r.esuppidx;
      %suppidx(ii,bm,2)=r.suppidx2;
      %esuppidx(ii,bm,2)=r.esuppidx2;
      
      % separate spatial, or, sf correlations -- too messy?
      tunecorr(bm,bm,ii,1)=r.spc;
      etunecorr(bm,bm,ii,1)=r.espc;
      tunecorr(bm,bm,ii,2)=r.orc;
      etunecorr(bm,bm,ii,2)=r.eorc;
      tunecorr(bm,bm,ii,3)=r.sfc;
      etunecorr(bm,bm,ii,3)=r.esfc;
      
      suppidx(:,bm,ii)=r.suppidx;
      esuppidx(:,bm,ii)=r.esuppidx;
      tempresp(:,bm,ii)=r.tempresp;
   end
end

expstr={'Rev','NR','GR'};
instr={'pix','pow','psf'};
cnfstr=expstr;
outstr={'lin','thr','ethr','full'};

expidx1=[1 1 2 ];
expidx2=[3 2 3 ];
stidx1=[4 4 2 ];
stidx2=[2 3 3 ];
cnfidx1=[1 1 3 ];
cnfidx2=[2 3 2 ];
nloutidx1=[2 2 2 ];

sunits='xc';
axisrange=[-0.2 1 0 15];

DOSPACE=2;

rowcount=length(expidx1);
colcount=3;

SEC=2.0;

fprintf('DOSPACE=%d SEC=%.1f\n',DOSPACE,SEC);

figure(1);

% figure out some stuff about how to group cells
goodidx4=find(~isnan(predxc(:,1)) & ~isnan(predxc(:,2)) & ...
              ~isnan(predxc(:,3)));
%goodidx2=goodidx3;


cumspacecomp=zeros(length(goodidx4),3,rowcount);


for ii=1:rowcount,
   
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
      set1=squeeze(spacecorr2(expidx1(ii),expidx2(ii),:,:));
      set2=squeeze(spacecorr2(expidx1(ii),expidx2(ii),:,:)+...
                   espacecorr2(expidx1(ii),expidx2(ii),:,:).*SEC);
      nroot=sprintf('%s/%s space2',...
                    expstr{expidx1(ii)},expstr{expidx2(ii)});
   elseif DOSPACE==3,
      % pull out temp corr data
      set1=squeeze(spacecorr3(expidx1(ii),expidx2(ii),:,:));
      set2=squeeze(spacecorr3(expidx1(ii),expidx2(ii),:,:)+...
                   espacecorr3(expidx1(ii),expidx2(ii),:,:).*SEC);
      nroot=sprintf('%s/%s space3',...
                    expstr{expidx1(ii)},expstr{expidx2(ii)});
   elseif DOSPACE==4,
      % pull out full/or/sf corr data
      set1=squeeze(tunecorr(expidx1(ii),expidx2(ii),:,:));
      set2=squeeze(tunecorr(expidx1(ii),expidx2(ii),:,:)+...
                   etunecorr(expidx1(ii),expidx2(ii),:,:).*SEC);
      nroot=sprintf('%s/%s or',...
                    expstr{expidx1(ii)},expstr{expidx2(ii)});
   else
      
   end
   nastr={[nroot,' all sig'],[nroot,' pos sig'],[nroot,' neg sig']};
   nbstr={'ns','ns','ns'};
   
   set1(find(set1>0.999))=0.999;
   cumspacecomp(:,:,ii)=set1(goodidx4,:);
   
   % figure out which cells have appropriate est data
   p1=predxc(:,expidx1(ii),1,cnfidx1(ii),nloutidx1(ii));
   p2=predxc(:,expidx2(ii),1,cnfidx2(ii),nloutidx1(ii));
   goodidx=find(~isnan(p1) & ~isnan(p2) & ~isnan(set1(:,1)));
   goodidx=intersect(goodidx,goodidx2);
   p1=p1(goodidx);
   p2=p2(goodidx);
   set1=set1(goodidx,:);
   set2=set2(goodidx,:);
   
   % classify as all-three or not
   goodidx3=find(sum(isnan(predxc(goodidx,:)),2)==0);
   d1=set1(goodidx3,:);
   d2=set1(setdiff(1:size(set1,1),goodidx3),:);
   
   
   % histogram full, pos, negative component correlations between
   % strfs of the two classes of interest
   for jj=1:3,

      % classify as to whether spatial correlations are
      % significantly different
      d1=set1(find(set2(:,jj)<1),jj);
      d2=set1(find(set2(:,jj)>=1),jj);
      
      % classify as to review st-sep kernel is better than hybrid
      p1=predxcst(goodidx,1,1,stidx1(ii));
      p2=predxcst(goodidx,1,1,stidx2(ii));
      e1=prederrst(goodidx,1,1,stidx1(ii));
      e2=prederrst(goodidx,1,1,stidx2(ii));
      
      d1=set1(find(p1-e1-e2 >= p2),jj); %-e1-e2
      d2=set1(find(p1-e1-e2 < p2),jj);
      
      % classify as either sig or ns based on pred accuracy
      e1=prederr(goodidx,expidx1(ii),1,cnfidx1(ii),nloutidx1(ii));
      e2=prederr(goodidx,expidx2(ii),1,cnfidx2(ii),nloutidx1(ii));
      d1=set1(find(p1-e1.*2>=0 & p2-e2.*2>=0),jj);
      d2=set1(find(p1-e1.*2<0 | p2-e2.*2<0),jj);
      
      
      subplot(rowcount,colcount,jj+(ii-1)*colcount);
      [p1,p2]=histcomp(d1(:),d2(:),nastr{jj},nbstr{jj},...
                       sunits,axisrange);
   end
   %subplot(rowcount,colcount,ii*colcount);
   %plotcomp([d1(:,2);d2(:,2)],[d1(:,3);d2(:,3)],nastr{2},nastr{3},...
   %         [-.2 1 -.2 1],[ones(size(d1,1),1);2*ones(size(d2,1),1)]);
end

TREP=500;
TTAIL=0;
[pa,ma]=randpairtest(cumspacecomp(:,1,1),cumspacecomp(:,1,2),TREP,TTAIL);
[pp,mp]=randpairtest(cumspacecomp(:,2,1),cumspacecomp(:,2,2),TREP,TTAIL);
[pn,mn]=randpairtest(cumspacecomp(:,3,1),cumspacecomp(:,3,2),TREP,TTAIL);

%[pa ma; pp mp; pn mn]

colormap(gray);
set(gcf,'PaperOrientation','landscape','PaperPosition',[0.25 0.25 10.5 8]);

figure(2);
clf
for ii=1:rowcount,
   
   subplot(rowcount,3,ii*3-2);
   
   set1=squeeze(tempresp(:,ii,:));
   goodidx=find(~isnan(set1(1,:)));
   %goodidx=intersect(goodidx,goodidx3);
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
   
   title(sprintf('%s mean tempresp (n=%d)',expstr{ii},size(set1,2)));
   if ii==rowcount,
      xlabel('time lag (ms)');
   end

   subplot(rowcount,3,ii*3-1);
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
   catidx(find(p1-e1.*2<0 | p2-e2.*2<0))=2;
   
   plotcomp(set1,set2,[expstr{expidx1(ii)},' adapt'],expstr{expidx2(ii)},...
            [0 1 0 1],catidx);
end

colormap(gray);
set(gcf,'PaperOrientation','landscape','PaperPosition',[0.25 0.25 10.5 8]);
 
figure(3);
clf
for ii=1:rowcount,
   
   subplot(rowcount,3,ii*3-2);
   
   % mean suppidx over time
   set1=squeeze(suppidx(3:end,ii,:));
   goodidx=find(~isnan(set1(1,:)));
   %goodidx=intersect(goodidx,goodidx3);
   set1=set1(:,goodidx);
   
   plot([7 (size(set1,1)*14-7)]+28,[0 0],'k--');
   hold on
   plot((7:14:(size(set1,1)*14-7))+28,mean(set1,2));
   hold off
   axis([0 size(set1,1)*14+28 0.2 0.6]);
   title(sprintf('%s mean suppidx (n=%d)',expstr{ii},size(set1,2)));
   if ii==rowcount,
      xlabel('time lag (ms)');
   end
   
   % suppidx at one point compared across classes
   subplot(rowcount,3,ii*3-1);
   set1=squeeze(suppidx(1,expidx1(ii),:));
   set2=squeeze(suppidx(1,expidx2(ii),:));
   goodidx=find(~isnan(set1) & ~isnan(set2));
   goodidx=intersect(goodidx,goodidx2);
   set1=set1(goodidx);
   set2=set2(goodidx);
   
   p1=predxc(goodidx,expidx1(ii),1,cnfidx1(ii),nloutidx1(ii));
   p2=predxc(goodidx,expidx2(ii),1,cnfidx2(ii),nloutidx1(ii));
   e1=prederr(goodidx,expidx1(ii),1,cnfidx1(ii),nloutidx1(ii));
   e2=prederr(goodidx,expidx2(ii),1,cnfidx2(ii),nloutidx1(ii));
   
   catidx=ones(length(goodidx),1).*1;
   catidx(find(p1-e1.*2<0 | p2-e2.*2<0))=2;
   
   plotcomp(set1,set2,[expstr{expidx1(ii)},' supp2'],expstr{expidx2(ii)},...
            [0 1 0 1],catidx);
   
   % suppidx at one point compared across classes
   subplot(rowcount,3,ii*3);
   set1=squeeze(suppidx(2,expidx1(ii),:));
   set2=squeeze(suppidx(2,expidx2(ii),:));
   goodidx=find(~isnan(set1) & ~isnan(set2));
   goodidx=intersect(goodidx,goodidx2);
   set1=set1(goodidx);
   set2=set2(goodidx);
   
   p1=predxc(goodidx,expidx1(ii),1,cnfidx1(ii),nloutidx1(ii));
   p2=predxc(goodidx,expidx2(ii),1,cnfidx2(ii),nloutidx1(ii));
   e1=prederr(goodidx,expidx1(ii),1,cnfidx1(ii),nloutidx1(ii));
   e2=prederr(goodidx,expidx2(ii),1,cnfidx2(ii),nloutidx1(ii));
   
   catidx=ones(length(goodidx),1).*1;
   catidx(find(p1-e1.*2<0 | p2-e2.*2<0))=2;
   
   plotcomp(set1,set2,[expstr{expidx1(ii)},' supp2'],expstr{expidx2(ii)},...
            [0 1 0 1],catidx);
end

colormap(gray);
set(gcf,'PaperOrientation','landscape','PaperPosition',[0.25 0.25 10.5 8]);
 

return

keyboard



MINLEN=500;
PMIN=0.01;
expstr={'Rev','GR','NR'};
instr={'pix','pow','psf'};
cnfstr=expstr;
outstr={'lin','thr','ethr','full'};





expidx1=[1 1 ];
expidx2=[2 3 ];
cnfidx=[1 1 ];
nloutidx1=[2 2 ];
nloutidx2=[2 2 ];

for ii=1:length(expidx1(:)),
   subplot(rowcount,colcount,ii);
   n1=[expstr{expidx1(ii)},'-',cnfstr{cnfidx(ii)},'-',instr{nlinidx1(ii)},...
       '-',outstr{nloutidx1(ii)}];
   n2=[expstr{expidx2(ii)},'-',cnfstr{cnfidx(ii)},'-',instr{nlinidx2(ii)},...
       '-',outstr{nloutidx2(ii)}];
   set1=predxc(:,expidx1(ii),nlinidx1(ii),cnfidx(ii),nloutidx1(ii));
   set2=predxc(:,expidx2(ii),nlinidx2(ii),cnfidx(ii),nloutidx2(ii));
   goodidx=find(~isnan(set1) & ~isnan(set2));
   set1=set1(goodidx);
   set2=set2(goodidx);
   p1=predp(goodidx,expidx1(ii),nlinidx1(ii),cnfidx(ii),nloutidx1(ii));
   p2=predp(goodidx,expidx2(ii),nlinidx2(ii),cnfidx(ii),nloutidx2(ii));
   e1=prederr(goodidx,expidx1(ii),nlinidx1(ii),cnfidx(ii),nloutidx1(ii));
   e2=prederr(goodidx,expidx2(ii),nlinidx2(ii),cnfidx(ii),nloutidx2(ii));
   
   set1hi=find(set1-e1 > set2+e2);
   set2hi=find(set1+e1 < set2-e2);
   sigidx=ones(size(set1))*2;
   sigidx(set1hi)=3;
   sigidx(set2hi)=1;
   
   %sigidx=2-(set1(goodidx)+e1(goodidx) < set2(goodidx)-e2(goodidx));
   %sigidx=2-(p1<PMIN | p2<PMIN);
   
   [p,m]=plotcomp(set1,set2,n1,n2,[-0.4 1.0 -0.4 1.0],sigidx);

   fprintf('%s v %s: %.2f v %.2f (p<%.3f) %d > %d < %d / %d \n',...
           n1,n2,mean(set1),mean(set2),p(end),sum(sigidx==3),...
           sum(sigidx==2),sum(sigidx==1),length(sigidx))
   fprintf('                    %.2f v %.2f\n',...
           mean(set1.*abs(set1)),mean(set2.*abs(set2)))
end





