% function [rmax,alpha,beta]=singletrialceiling(r,[alpha],[beta],N[=100],p[=0.95]);
%
% calculate noise ceiling on correlation coefficient for
% single-trial spike data. to correct measured r^2 for single trial
% noise, use:
%
%             r2_corrected = r2_measured ./ (rmax.^2)
%
% inputs:
% r- vector of spikes
% alpha,beta - parameters of gamma function that produced r. if not
%              provided, these are fit by reversegamma
% N - number of random modulator functions to test
% p - rank of cc to return from among the random tests
%     [0.95]. smaller p will return larger rmax
%
% output:
% rmax - upper bound on correlation coefficient between single
%        trial data (provided) and theoretical modulator function 
%        (determined by gamma parameters)
% alpha,beta - gamma parameters (fit or just return inputs)
%
% created SVD 2004-04-26
%
function [rmax,alpha,beta]=singletrialceiling(r,alpha,beta,N,p);

% make sure r is all integers
r=round(r);

if ~exist('beta','var'),
   % fit gamma parameters for the single trial data
   [mu,alpha,beta]=reversegamma(r);
end

if ~exist('N','var'),
   N=100;
end
if ~exist('p','var'),
   p=0.95;
end

cc=zeros(N,1);
for ii=1:N,
   trp=reversegamma(r,alpha,beta);
   if var(r)>0 & var(trp)>0,
      cc(ii)=xcov(r,trp,0,'coeff');
   end
end

cc=sort(cc);
rmax=cc(round(p*N));


