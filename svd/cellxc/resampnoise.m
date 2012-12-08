% function rr=resampnoise(r,n)
%
% shuffled resampling
%
% assume r is 2D (time X space)... add n noise shuffles along third dim
%
function rr=resampnoise(r,n)

if not(exist('n')),
   n=1;
end
rlen=size(r,2);
respcount=size(r,2);
rr=zeros(rlen,respcount,n-1);
fprintf('Resampling response (length=%d, resamples=%d)\n',rlen,n);
for respidx=1:respcount,
   goodidx=find(~isnan(r(:,respidx)));
   for ii=1:n-1,
      tt=rand(length(goodidx),1);
      [stt,tidx]=sort(tt);
      rr(goodidx,respidx,ii)=r(goodidx(tidx));
   end
   
   badidx=find(isnan(r(:,respidx)));
   rr(badidx,respidx,:)=nan;
end

rr=cat(3,r,rr);





