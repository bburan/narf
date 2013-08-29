% function rsingle=xcnoiseceiling(pred,r_raster,N[=100]);
%
% Calculate noise ceiling on correlation coefficient for single-trial
% spike data, using methodology based on ideas of FET.  To correct
% measured r^2 for single trial noise, use:
%
%             xc_norm = xc_single ./ (ac_single.^2)
%
% inputs:
% pred - T x 1 vector of predicted firing rate
% r_raster - T x M matrix of M single-trial actual responses.
% N - number of random single-trial correlations to measure. 100
% should be fine in most cases
%
% output:
% xc_norm - upper bound on correlation coefficient between single
%           trial data for given p value
% ac_single - mean correlation between single trial data
%
% created SVD 2013-08-22
%
function [r_norm,rmax]=xcnoiseceiling(pred,r_raster,N,p);

M=size(r_raster,2);
validbins=find(~isnan(pred) & ~isnan(r_raster(:,1)));

mc=round(M./2);
sh=shuffle(1:M);
r1=r_raster(validbins,sh(1:mc));
r2=r_raster(validbins,sh((mc+1):end));

rM=xcov(pred(validbins),mean(r_raster(validbins,:),2),0,'coeff');
r12=xcov(mean(r1,2),mean(r2,2),0,'coeff');

rmax=sqrt(1./(1+0.5*(-M+M./r12)));

xc1=zeros(M,1);
for ii=1:M,
    if var(pred(validbins))>0 & var(r_raster(validbins,ii))>0,
        xc1(ii)=xcov(pred(validbins),r_raster(validbins,ii),0,'coeff');
    end
end
xc1mean=mean(xc1);
r_norm=xc1mean./rmax;


return



if ~exist('N','var'),
   N=100;
end
if ~exist('p','var'),
   p=0.95;
end

% make sure r is all integers
r_raster=round(r_raster);

r=r_raster(:);  % make a vector

[mu,alpha,beta]=reversepoisson(r);

if 1
    rmax=singletrialceiling(r,alpha,beta);
else
    cc=zeros(N,1);
    for ii=1:N,
        trp=reversegamma(r,alpha,beta);
        if var(r)>0 & var(trp)>0,
            cc(ii)=xcov(r,trp,0,'coeff');
        end
    end
    cc=sort(cc);
    rmax=cc(round(p*N));
end

M=size(r_raster,2);
validbins=find(~isnan(pred) & ~isnan(r_raster(:,1)));
xc1=zeros(M,1);
for ii=1:M,
    if var(pred(validbins))>0 & var(r_raster(validbins,ii))>0,
        xc1(ii)=xcov(pred(validbins),r_raster(validbins,ii),0,'coeff');
    end
end

r_norm=sqrt(mean(xc1.^2)./rmax.^2);
return
r_single=mean(xc1);
r_norm=r_single./rmax;

%keyboard

