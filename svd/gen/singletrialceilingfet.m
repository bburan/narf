% function rsingle=singletrialceilingfet(r,N[=100],p[=0.95]);
%
% calculate noise ceiling on correlation coefficient for
% single-trial spike data. to correct measured r^2 for single trial
% noise, use:
%
%             r2_corrected = r2_measured ./ (rmax.^2)
%
% inputs:
% r- vector of spikes
% N - number of random modulator functions to test
% p - rank of cc to return from among the random tests
%     [0.95]. smaller p will return larger rmax
%
% output:
% rmax - upper bound on correlation coefficient between single
%        trial data (provided) and theoretical modulator function 
%        (determined by gamma parameters)
%
% created SVD 2012-02-06
%
function rsingle=singletrialceilingfet(r,N,p);

if ~exist('N','var'),
   N=100;
end
if ~exist('p','var'),
   p=0.95;
end

cc=zeros(N,1);
repcount=size(r,2);
if repcount<=1,
    error('more than one rep required');
end
for ii=1:N,
    x=ceil(rand*repcount);
    y=ceil(rand*(repcount-1));
    if y>=x,
        y=y+1;
    end
    if var(r(:,x))>0 & var(r(:,y))>0,
      cc(ii)=xcov(r(:,x),r(:,y),0,'coeff');
   end
end

rsingle=mean(cc);
return

cc=sort(cc);
rmax=cc(round(p*N));


