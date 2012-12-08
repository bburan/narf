% function rr=respresamptime(r,n,skipfirstcol,usecolavg)
%
% bootstrap resampling
%
% resample binned spike function, r, n times, excluding 1/n of the
% data in each sample.  for respresampTIME, excluded segment is always
% a segement of the movie, rather than trials (ie, for multi-trial
% data sets with more than one column in r)
%
function rr=respresamptime(r,n,skipfirstcol,usecolavg)

if length(find(isnan(r)))>0,
   usenan=1;
else
   usenan=0;
end

if not(exist('n')),
   n=1;
end
if not(exist('skipfirstcol')),
   skipfirstcol=0;
end
if not(exist('usecolavg')),
   usecolavg=1;
end
if skipfirstcol,
   n=n-1;
end

rlen=size(r,1);
rcount=size(r,2);
disp(sprintf('Resampling response (length=%d, count=%d, resamples=%d)',...
             rlen,rcount,n));
   

if usecolavg,
   if usenan,
      rr=ones(rlen,n)*nan;
   else
      rr=-ones(rlen,n);
   end   
   
   % dispose of separate trials
   if rcount>1,
      % find mean across trials
      rtemp=r;
      rnz=rcount-sum(rtemp==-1,2);
      rtemp=rtemp.*(rtemp>0);
      
      rs=sum(rtemp,2)./(rnz+(rnz==0));
      rs(find(rnz==0))=-1;
   else
      rs=r;
   end
   
   % break valid frames into n equal length segments
   rok=find(rs>=0);
   rendidx=rok(round((1:n)' .* length(rok)./n));
   rstartidx=[1;rendidx(1:n-1)+1];
   
   % resampled response is the mean response minus one of the segments
   for ii=1:n,
      rr(1:rstartidx(ii)-1,ii)=rs(1:rstartidx(ii)-1);
      rr(rendidx(ii)+1:rlen,ii)=rs(rendidx(ii)+1:rlen);
   end
   
   if skipfirstcol,
      rr=[rs rr];
   end
else
   rr=-ones(rlen,rcount,n);
   
   rok=sum(~isnan(r),2)>0;
   rendidx=rok(round((1:n)' .* length(rok)./n));
   
end
