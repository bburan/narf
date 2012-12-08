% function crfcomp(cellid,batchid);
%
function crfcomp(cellid,batchid);

cellfiledata=cellfiletimes(cellid,batchid);
stimfile=[cellfiledata(1).path,cellfiledata(1).stimfile];
respfile=[cellfiledata(1).path,cellfiledata(1).respfile];

stim=loadimfile(stimfile,1,0,16,16,16);
resp=resploadatt(respfile,'pfth',1,0,1);
resp=respvarbins(resp,125,275);

mstim=mean(stim(:));
m=squeeze(mean(mean(stim,1),2));

c=stim.^2;
c=sqrt(squeeze(mean(mean(c,1),2))-m.^2);

% remove fixations that aren't associated with a specific
% attention condition
respidx=1;
SKIPATT1=1;

if SKIPATT1,
   okidx=find(sum(~isnan(resp(:,1,2:end)),3)>0);
   resp=resp(okidx,:,2:end);
   c=c(okidx);
   m=m(okidx);
   %resp=resp*1000; % convert to Hz
end

rlen=size(resp,1);
respcount=size(resp,2);
attcount=size(resp,3);
bincount=8;

h=zeros(bincount,attcount);
mmin=0;
mmax=nanmax(resp(:))-nanstd(resp(:));
edges=linspace(mmin,mmax,bincount);
edges(end)=inf;
for attidx=1:attcount,
   binuse=find(~isnan(resp(:,respidx,attidx)));
   h(:,attidx)=histc(resp(binuse,respidx,attidx),edges);
   h(:,attidx)=h(:,attidx)./sum(h(:,attidx));
end

%compute crf

[sc,sic]=sort(c);
[sm,sim]=sort(m);
blen=length(c);
medges=zeros(bincount+1,1);      
cedges=zeros(bincount+1,1);      
mmid=zeros(bincount,1);      
cmid=zeros(bincount,1);
for bb=1:bincount,
   % [ round((bb-1)*blen/bincount+1) round(bb*blen/bincount+1) ]
   
   medges(bb)=m(sim(round((bb-1)*blen/bincount+1)));
   mmid(bb)=mean(m(sim(round((bb-1)*blen/bincount+1):...
                       round(bb*blen/bincount))));
   cedges(bb)=c(sic(round((bb-1)*blen/bincount+1)));
   cmid(bb)=mean(c(sic(round((bb-1)*blen/bincount+1):...
                       round(bb*blen/bincount))));
end
medges(end)=inf;
cedges(end)=inf;

hc=zeros(bincount,attcount);
hce=zeros(bincount,attcount);
gc=zeros(bincount,attcount);
hm=zeros(bincount,attcount);
hme=zeros(bincount,attcount);
csig=zeros(4,attcount);
msig=zeros(4,attcount);
for attidx=1:attcount,
   for binidx=1:bincount,
      cbins=find(c>=cedges(binidx) & c<cedges(binidx+1) & ...
                 ~isnan(resp(:,respidx,attidx)));
      
      hc(binidx,attidx)=mean(resp(cbins,respidx,attidx));
      hce(binidx,attidx)=std(resp(cbins,respidx,attidx))./sqrt(length(cbins));
      
      %tr=resp(cbins,respidx,attidx)-mean(resp(cbins,respidx,attidx));
      %ts=c(cbins)-mean(c(cbins));
      tr=resp(cbins,respidx,attidx)-nanmean(resp(:,respidx,attidx));
      ts=c(cbins)-mean(c(:));
      %gc(binidx,attidx)=mean(tr.*ts)./mean(ts.^2);
      
      
      mbins=find(m>=medges(binidx) & m<medges(binidx+1) & ...
                 ~isnan(resp(:,respidx,attidx)));
      
      hm(binidx,attidx)=mean(resp(mbins,respidx,attidx));
      hme(binidx,attidx)=std(resp(mbins,respidx,attidx))./sqrt(length(mbins));
      
   end
   
   ridx=find(~isnan(resp(:,respidx,attidx)));
   
   csig(:,attidx)=fitsigmoid(c(ridx),resp(ridx,respidx,attidx),0)';
   %msig(:,attidx)=fitsigmoid(m(ridx),resp(ridx,respidx,attidx),0)';
   
end


sfmt={'bo','ro','go','ko'};
ffmt={'b-','r-','g-','k-'};

figure(1);
clf

subplot(3,1,1);
redges=edges';
redges(end)=redges(end-1)*2-redges(end-2);
bar(redges,h);
%hold on
%plot(redges,flipud(cumsum(flipud(h),1)));
%hold off
title([cellid,' response histogram']);
legend('a1','a2','a3','a4');

subplot(3,1,2);
for attidx=1:attcount
   errorbar(cmid,hc(:,attidx),hce(:,attidx),sfmt{attidx});
   hold on
   plot(cmid,sigmoid(csig(:,attidx),cmid),ffmt{attidx});
end
hold off
title([cellid,' contrast gain']);

ttext=sprintf('(%.1f,%.2f,%.3f,%.3f)\n',csig);
ttext=sprintf('(smean,slope,rmin,rmax)\n%s',ttext);
a=axis;
x1=a(1)+(a(2)-a(1))./100;
y1=a(4)-(a(4)-a(3))./100;

text(x1,y1,ttext,'VerticalAlignment','top');


subplot(3,1,3);
hl=zeros(attcount,1);
for attidx=1:attcount
   ll=errorbar(mmid,hm(:,attidx),hme(:,attidx),sfmt{attidx});
   hl(attidx)=ll(2);
   hold on
   plot(mmid,hm(:,attidx),ffmt{attidx});
end
hold off
title('mean lum vs response');

legend(hl,'a1','a2','a3','a4');
set(gcf,'PaperOrientation','Portrait','PaperPosition',[0.25 0.25 8 10.5]);

%keyboard


return

sql='SELECT * FROM sRunData where batch=11 ORDER BY cellid;';

dbopen;
rundata=mysql(sql);

for ii=3:length(rundata),
   crfcomp(rundata(ii).cellid,11);
   drawnow;
   print -f1 -Pgcolor
end
