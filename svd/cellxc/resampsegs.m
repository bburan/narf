% function [rstartidx,rendidx]=resampsegs(r,n)
%
% bootstrap resampling
%
% resample binned spike function, r, n times, excluding 1/n of the
% data in each sample.  for respresampTIME, excluded segment is always
% a segement of the movie, rather than trials (ie, for multi-trial
% data sets with more than one column in r)
%
function [rstartidx,rendidx]=resampsegs(r,n)

if not(exist('n')),
   n=1;
end

respcount=size(r,2);
rstartidx=zeros(n,respcount);
rendidx=zeros(n,respcount);

fprintf('Resampling response (count=%d): ',n);
for respidx=1:respcount,
   rr=r(:,respidx);   % only look at first column of r
   rlen=length(rr);
   
   % figure out valid time bins
   rok=find(~isnan(rr));
   
   if length(rok)>=n,
      fprintf(' (r=%d len=%d)',respidx,length(rok));
      
      % figure out start and stop times of n equal length segments
      rendidx(:,respidx)=rok(round((1:n)' .* length(rok)./n));
      rstartidx(:,respidx)=[1;rendidx(1:n-1,respidx)+1];
      
   elseif length(rok)==0,
     fprintf(' (%d len=0!!)',respidx);
     rendidx(:,respidx)=zeros(n,1);
     rstartidx(:,respidx)=ones(n,1);
     %[rstartidx rendidx]
   else
     fprintf(' (%d len=%d < %d resamps)',respidx,length(rok),n);
     rendidx(:,respidx)=[rok;ones(n-length(rok),1).*(rok(end)-1)];
     rstartidx(:,respidx)=[rok;ones(n-length(rok),1).*rok(end)];
     %[rstartidx rendidx]
   end
      
end
fprintf('\n');



