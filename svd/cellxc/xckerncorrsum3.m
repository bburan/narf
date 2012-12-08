
batches=[23; 29; 32];
batches2=[86; 83; 85; 84; 88];

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

figure(1);
clf

for ii=1:compcount,
   classexists=~isnan(predxc(:,expidx1(ii))) & ~isnan(predxc(:,expidx2(ii)));
   allclassexists=sum(~isnan(predxc(:,:,1)),2)==3;
   
   p1=predxc(:,expidx1(ii),1,cnfidx1(ii),nloutidx1(ii));
   p2=predxc(:,expidx2(ii),1,cnfidx2(ii),nloutidx1(ii));
   e1=prederr(:,expidx1(ii),1,cnfidx1(ii),nloutidx1(ii));
   e2=prederr(:,expidx2(ii),1,cnfidx2(ii),nloutidx1(ii));
   
   %withingoodidx=find(classexists & p1-e1.*2>0 & p2-e2.*2>0);
   withingoodidx=find(classexists);
   allwithingoodidx=find(classexists & p1-e1.*2>0 & p2-e2.*2>0);
   %allwithingoodidx=find(allclassexists & p1-e1.*2>0 & p2-e2.*2>0);
   %allwithingoodidx=find(allclassexists);
   
   pr1=predxc(:,expidx1(ii),1,1,nloutidx1(ii));
   pr2=predxc(:,expidx2(ii),1,1,nloutidx1(ii));
   pst1=predxcst(:,1,1,stcnfidx1(ii),stnlidx1(ii));
   pst2=predxcst(:,1,1,stcnfidx2(ii),stnlidx2(ii));
   ppos1=predxcst(:,2,1,stcnfidx1(ii),posnlidx1(ii));
   ppos2=predxcst(:,2,1,stcnfidx2(ii),posnlidx2(ii));
   pneg1=predxcst(:,3,1,stcnfidx1(ii),posnlidx1(ii));
   pneg2=predxcst(:,3,1,stcnfidx2(ii),posnlidx2(ii));

   
   pall=[pr1 pst1 ppos1 pneg1 pneg2 ppos2 pst2 pr2];
   
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


