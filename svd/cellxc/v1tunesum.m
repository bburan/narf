% v1tunesum.m
%
% summarize tuning over jackknifed kernels for different stim
% classes
%
% created SVD 10/2/03
%

clear

dbopen;

figure(1);
figure(2);

SPACE=0;

if SPACE,
   batches=[52; ...   % rev: 1.0 0.95 0.9 0.85
            56; ...   % gr:  ""
            60];      % nr:  ""
elseif 1,
   batches=[51; ...   % rev: 1.0 0.95 0.9 0.85
            55; ...   % gr:  ""
            59];      % nr:  ""
end

yn=input('reload and re-calc tuning (y/[n])? ','s');
fname=sprintf('/auto/k5/david/tmp/v1tunsum.s%d.mat',SPACE);

if (length(yn)>0 & yn(1)=='y') | ~exist(fname,'file'),
   
   s=sprintf('%d,',batches);
   s(end)=')';
   
   sql=['SELECT cellid,count(id) as batchcount FROM sRunData',...
        ' WHERE batch in (',s,' GROUP BY cellid ORDER BY cellid'];
   resdata=mysql(sql);
   
   clear tunedata fitdata
   
   cellcount=length(resdata);
   or=zeros(cellcount,length(batches)).*nan;
   orw=zeros(cellcount,length(batches)).*nan;
   sf=zeros(cellcount,length(batches)).*nan;
   sfw=zeros(cellcount,length(batches)).*nan;
   lat=zeros(cellcount,length(batches)).*nan;
   arat=zeros(cellcount,length(batches)).*nan;
   eor=zeros(cellcount,length(batches)).*nan;
   eorw=zeros(cellcount,length(batches)).*nan;
   esf=zeros(cellcount,length(batches)).*nan;
   esfw=zeros(cellcount,length(batches)).*nan;
   elat=zeros(cellcount,length(batches)).*nan;
   earat=zeros(cellcount,length(batches)).*nan;
   
   for cellidx=1:length(resdata),
      [tunedata,tfitdata]=xcshowtune(resdata(cellidx).cellid,batches);
      for ii=1:length(tfitdata),
         bid=find(tfitdata(ii).batch==batches);
         fitdata(cellidx,bid)=tfitdata(ii);
         
         or(cellidx,bid)=tfitdata(ii).orfit(1);
         orw(cellidx,bid)=tfitdata(ii).orfit(2);
         sf(cellidx,bid)=tfitdata(ii).sffit(1);
         sfw(cellidx,bid)=tfitdata(ii).sffit(2);
         lat(cellidx,bid)=tfitdata(ii).timefit(1);
         arat(cellidx,bid)=1-tfitdata(ii).timefit(6)./...
             tfitdata(ii).timefit(3);
         eor(cellidx,bid)=tfitdata(ii).orfiterr(1);
         eorw(cellidx,bid)=tfitdata(ii).orfiterr(2);
         esf(cellidx,bid)=tfitdata(ii).sffiterr(1);
         esfw(cellidx,bid)=tfitdata(ii).sffiterr(2);
         elat(cellidx,bid)=tfitdata(ii).timefiterr(1);
         earat(cellidx,bid)=1-tfitdata(ii).timefiterr(6)./...
             tfitdata(ii).timefiterr(3);
         
      end
      close all
   end
   
   save(fname);
else
   load(fname);
end

batchcount=length(batches);
paircount=batchcount*(batchcount-1)./2;

p1=[1 1 2 1 2 3];
p2=[2 3 3 4 4 4];

PLOTERR=1;

figure(1);
for pidx=1:paircount,
   subplot(6,paircount,pidx);
   if PLOTERR,
      errscat(or(:,p1(pidx)),or(:,p2(pidx)),...
              eor(:,p1(pidx)),eor(:,p2(pidx)));
   else
      scatter(or(:,p1(pidx)),or(:,p2(pidx)),'.');
   end
   
   hold on
   plot([0 180],[0 180],'--');
   hold off
   title(sprintf('or: %d v %d',...
                 batches(p1(pidx)),batches(p2(pidx))));
   axis([0 180 0 180]);
   axis square
   
   subplot(6,paircount,pidx+paircount);
   if PLOTERR,
      errscat(orw(:,p1(pidx)),orw(:,p2(pidx)),...
              eorw(:,p1(pidx)),eorw(:,p2(pidx)));
   else
      scatter(orw(:,p1(pidx)),orw(:,p2(pidx)),'.');
   end
   hold on
   plot([0 60],[0 60],'--');
   hold off
   title(sprintf('orw: %d v %d',...
                 batches(p1(pidx)),batches(p2(pidx))));
   axis([0 100 0 100]);
   axis square
   
   subplot(6,paircount,pidx+paircount.*2);
   if PLOTERR,
      errscat(sf(:,p1(pidx)),sf(:,p2(pidx)),...
              esf(:,p1(pidx)),esf(:,p2(pidx)));
   else
      scatter(sf(:,p1(pidx)),sf(:,p2(pidx)),'.');
   end
   hold on
   plot([0 4],[0 4],'--');
   hold off
   title(sprintf('sf: %d v %d',...
                 batches(p1(pidx)),batches(p2(pidx))));
   axis([0 4 0 4]);
   axis square
   
   subplot(6,paircount,pidx+paircount.*3);
   if PLOTERR,
      errscat(sfw(:,p1(pidx)),sfw(:,p2(pidx)),...
              esfw(:,p1(pidx)),esfw(:,p2(pidx)));
   else
      scatter(sfw(:,p1(pidx)),sfw(:,p2(pidx)),'.');
   end
   hold on
   plot([0 3],[0 3],'--');
   hold off
   title(sprintf('sfw: %d v %d',...
                 batches(p1(pidx)),batches(p2(pidx))));
   axis([0 3 0 3]);
   axis square
   
   subplot(6,paircount,pidx+paircount.*4);
   if PLOTERR,
      errscat(lat(:,p1(pidx)),lat(:,p2(pidx)),...
              elat(:,p1(pidx)),elat(:,p2(pidx)));
   else
      scatter(lat(:,p1(pidx)),lat(:,p2(pidx)),'.');
   end
   hold on
   plot([0 100],[0 100],'--');
   hold off
   title(sprintf('lat: %d v %d',...
                 batches(p1(pidx)),batches(p2(pidx))));
   axis([0 100 0 100]);
   axis square
   
   subplot(6,paircount,pidx+paircount.*5);
   if PLOTERR,
      errscat(arat(:,p1(pidx)),arat(:,p2(pidx)),...
              earat(:,p1(pidx)),earat(:,p2(pidx)));
   else
      scatter(arat(:,p1(pidx)),arat(:,p2(pidx)),'.');
   end
   hold on
   plot([-1 1],[-1 1],'--');
   hold off
   title(sprintf('arat: %d v %d',...
                 batches(p1(pidx)),batches(p2(pidx))));
   axis([-1 1 -1 1]);
   axis square
   
end
