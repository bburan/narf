% function [m,se]=jackmse(x,y,n);
%
% computes MSE between each column of x and y.  (x and y should be
% matrices of the same size.)  use n jackknifes (default 20)
%
% m - mean
% se - jackknife estimate of standard error
%
function [m,se]=jackmse(x,y,n);

if ~exist('n','var'),
   n=20;
end

if length(x)==1,
   warning('jackmeanerr: length(x)=1, setting m=se=x');
   m=x;
   se=x;
   return
end

if n>size(x,1),
   n=size(x,1);
end

jackstep=size(x,1)/n;

mi=zeros(n,size(x,2));
for ii=1:n,
   jackrange=[1:round((ii-1)*jackstep) round(jackstep*ii+1):size(x,1)];
   mi(ii,:)=sqrt(mean((x(jackrange,:)-y(jackrange,:)).^2))./...
            std(y(jackrange,:));
end

m=nanmean(mi,1);
se=nanstd(mi,0) * sqrt(n-1);



