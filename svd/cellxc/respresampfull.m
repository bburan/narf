% function rr=respresampfull(r,n,includefull,usecolavg)
%
% bootstrap resampling
%
% resample binned spike function, r, n times, excluding 1/n of the
% data in each sample.  for respresampTIME, excluded segment is always
% a segement of the movie, rather than trials (ie, for multi-trial
% data sets with more than one column in r)
%
function rr=respresampfull(r,N,includefull,usecolavg)

if length(find(isnan(r)))>0,
   usenan=1;
else
   usenan=0;
   r(r==-1)=nan;
end

if not(exist('N')),
   N=1;
end
if not(exist('includefull')),
   includefull=0;
end
if not(exist('usecolavg')),
   usecolavg=1;
end
if includefull,
   n=N-1;
else
   n=N;
end

rlen=size(r,1);
rcount=size(r,2);
disp(sprintf('Resampling response (length=%d, count=%d, resamples=%d)',...
             rlen,rcount,n));

if usecolavg,
   % assume r is time X trial matrix
   rlen=size(r,1);
   rcount=size(r,2);
   
   % dispose of separate trials
   if rcount>1,
      % find mean across trials
      rtemp=r.*(~isnan(r));
      rnz=rcount-sum(isnan(r),2);
      
      rs=sum(rtemp,2)./(rnz+(rnz==0));
      rs(find(rnz==0))=nan;
   else
      rs=r;
   end
   
else
   rlen=size(r);
   rs=reshape(r,rlen(1),prod(rlen(2:end)));
end

rsize=size(rs);
rok=find(sum(~isnan(rs),2)>0);
   
if usenan,
   rr=ones([rsize n])*nan;
else
   rr=-ones([rsize,n]);
end

rendidx=rok(round((1:n)' .* length(rok)./n));
rstartidx=[1;rendidx(1:n-1)+1];

% resampled response is the mean response minus one of the segments
for ii=1:n,
   rr(1:rstartidx(ii)-1,:,ii)=rs(1:rstartidx(ii)-1,:);
   rr(rendidx(ii)+1:rlen,:,ii)=rs(rendidx(ii)+1:rlen,:);
end

if includefull,
   rr=cat(3,rs,rr);
end

rr=squeeze(reshape(rr,[rlen,N]));





