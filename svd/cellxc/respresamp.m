% function rr=respresamp(r,frac,n)
%
% resample binned spike function, r, n times, taking frac X100% of the
% spikes in each resampling.
%
function rr=respresamp(r,frac,n)

if not(exist('n')),
   n=length(frac);
end

rlen=size(r,1);
rcount=size(r,2);
rr=zeros(rlen,n);

fprintf('Resampling response (length=%d, count=%d, resamples=%d)\n',...
        rlen,rcount,n);

if rcount==1, % ie, single spike train, as in natrev/gratrev
   ri=find(r==1);
   % alternatively, ri=(1:rlen)';
   
   for ii=1:n,
      ruse=rand(rlen,1);
      if length(frac)==n,
         if frac(ii)==1,
	    rr(:,ii)=r;
         else
	    rr(:,ii)=r.*(ruse <= frac(ii));
         end
      else
         rr(:,ii)=r.*(ruse <= frac(1));
      end
   end
   
else  %rcount > 1, sample from different response sets
   % r=r.*(r>0);  % get rid of -1 terms
   usecount=ceil(frac.*(frac<1).*rcount);
   usetotal=sum(usecount);
   %useidx=mod((0:usetotal-1)',rcount)+1;
   
   % uniform distribution for picking random resamples
   useidx=ceil(rand(usetotal,1)*rcount);
   [y,ridx]=sort(rand(usetotal,1));
   jj=1;
   for ii=1:n,
      if frac(ii)==1,
         rtemp=r;
         rnz=rcount-sum(rtemp==-1,2);
      else
         rtemp=r(:,useidx(ridx(jj:jj+usecount(ii)-1)));
         jj=jj+usecount(ii);
         rnz=usecount(ii)-sum(rtemp==-1,2);
      end;
      rtemp=rtemp.*(rtemp>0);
      
      rr(:,ii)=sum(rtemp,2)./(rnz+(rnz==0));
      rr(find(rnz==0),ii)=-1;
   end
   
   %for ii=1:n,
   %   usecount=round(frac(ii)*rcount);
   %   [y,ridx]=sort(rand(rcount,1));
   %   rtemp=r(:,ridx(1:usecount));
   %   rnz=usecount-sum(rtemp==-1,2);
   %   rnz=rnz;
   %   rtemp=rtemp.*(rtemp>0);
   %   
   %   rr(:,ii)=sum(rtemp,2)./(rnz+(rnz==0));
   %   rr(find(rnz==0),ii)=-1;
   %end
end

