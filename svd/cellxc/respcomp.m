function respcomp(resp);


respidx=1;
SKIPATT1=1;

if SKIPATT1,
   okidx=find(sum(~isnan(resp(:,1,2:end)),3)>0);
   resp=resp(okidx,:,2:end);
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

figure(1);
clf
subplot(2,1,1);
plot(h);

subplot(2,1,2);
plot(flipud(cumsum(flipud(h),1)));

